"""Apply a perturbation file to a running 4DWCM build.

Whole_Cell_Minimal_Cell.py calls, only when --perturb is given:

    P = apply.Perturbation.from_file(path, head_directory)     # load, validate against the spec, resolve
    P.install(sim, sim_properties)                              # right after initSim, before anything reads the inputs
    P.after_particles(sim, sim_properties)                      # right after IC.initializeParticles
    P.report(sim, sim_properties)                               # before sim.run: read back what the simulation holds

Without --perturb none of this runs, so the unperturbed model is built exactly as before. The edits enter at the points
where the simulator takes the values in:

    promoter strengths         sim_properties['promoters']          (read by GIP.RNAP_binding and GIP.TranscriptionRate)
    initial protein / RNA      sim.distributeNumber                 (called by initializeProteins / MRNA / TRNA)
    initial metabolites, medium sim_properties['counts'] / ['medium']
    spreadsheet parameters     Rxns_ODE._EXCEL_CACHE                (every ODE / CME build reads the cached kinetic tables)
    RDME rate constants        sim.rateConst                        (shared, assembly and per-gene constants)
    CME transcription          GIP_rates.TranscriptionRate          (rebuilt every hook)
    diffusion                  sim.diffusionConst
    GIP formula constants      utility.GIP_rates module globals
    SMC floor                  sim_properties['smc_min']            (SpatialDnaDynamics._num_smc; 0 when P_0415 is knocked out)
"""

import json
import os
from collections import OrderedDict

import numpy as np

from . import perturbation as P
from .impact import gene_species

PER_GENE_SUFFIX = {'translation': '_translat', 'mrna_degradation': '_RNAdeg', 'secy_insertion': '_insertion', 'rnap_binding': 'RNAP_on'}
SMC = 'P_0415'


def _own_rna_dc(name):
    base = name[:-8] if name.endswith('_diffDna') else name[:-5] if name.endswith('_diff') else None
    return base is not None and base.startswith('R') and base != 'RNAP'


class Perturbation:
    def __init__(self, pert, spec, head, path=None):
        self.pert, self.spec, self.head, self.path = pert, spec, head.rstrip('/') + '/', path
        self.resolution = P.resolve(spec, pert)
        self.log = []                 # (what, requested, value handed to the simulator)
        self._placed = {}             # species -> total placed by the wrapped distributeNumber
        genes = spec['genes']
        self.ko = {k['gene']: k.get('mode', 'full') for k in pert.get('knockouts') or []}
        self.kd = {k['gene']: k['promoter_scale'] for k in pert.get('knockdowns') or []}
        # species whose t=0 placement changes: None = zero, int = new placed total
        self.place = {}
        for locus, mode in self.ko.items():
            if mode in ('full', 'initial_only'):
                for sid in gene_species(spec, locus):
                    self.place[sid] = 0
        for sid, v in (pert.get('initial_protein_counts') or {}).items():
            self.place[sid] = int(v)
        self.mrna_mean = dict(pert.get('initial_mrna_means') or {})
        # RDME rate constants by exact name, and per-gene scales by (template suffix, locus number or '*')
        self.rc = {e['constant']: e for e in pert.get('rdme_rate_constants') or []}
        self.gene_scale = {}
        for e in pert.get('gene_rate_scales') or []:
            num = '*' if e['gene'] == '*' else genes[e['gene']]['num']
            self.gene_scale[(e['rate'], num)] = e['scale']
        self.dc = {e['constant']: e['scale'] for e in pert.get('diffusion_scales') or []}

    # ------------------------------------------------------------------------------------------ construction
    @classmethod
    def from_file(cls, path, head):
        head = head.rstrip('/') + '/'
        pert = P.load(path)
        spec = load_spec(head)
        err, warn = P.validate(spec, pert)
        for w in warn:
            print('PERTURBATION warning: ' + w)
        if err:
            raise SystemExit('PERTURBATION %s is invalid:\n  ' % path + '\n  '.join(err))
        print('PERTURBATION %s: %d edits' % (pert.get('name'), len(P.resolve(spec, pert)['edits'])))
        return cls(pert, spec, head, path)

    def _rate_scale(self, rate, locus_num):
        s = 1.0
        for key in ((rate, '*'), (rate, locus_num)):
            if key in self.gene_scale:
                s *= self.gene_scale[key]
        return s

    # ------------------------------------------------------------------------------------------ install
    def install(self, sim, sim_properties):
        """Before any input is read: GIP constants, the cached kinetic tables, and wrappers on the sim."""
        import utility.GIP_rates as GIP
        from collections import OrderedDict as OD
        for e in self.pert.get('gip_constants') or []:
            self.log.append(('GIP_rates.%s' % e['constant'], getattr(GIP, e['constant']), e['value']))
            setattr(GIP, e['constant'], e['value'])
        if any(e['constant'] in ('ATPconc', 'UTPconc', 'CTPconc', 'GTPconc') for e in self.pert.get('gip_constants') or []):
            GIP.baseMap = OD({'A': GIP.ATPconc, 'U': GIP.UTPconc, 'G': GIP.GTPconc, 'C': GIP.CTPconc})   # built at import
        self._edit_tables()
        self._wrap_transcription(GIP)
        self._wrap_sim(sim)
        sim_properties['perturbation'] = self.pert.get('name')
        # the ODE code has the rate constants compiled in: build it in this run's own directory, never the shared tree
        import utility.Integrate as Integrate
        wd = sim_properties.get('working_directory')
        if wd:
            Integrate.ODE_BUILD_DIR = os.path.join(wd, 'ode_build')
        # the ODE forms of a fully knocked-out carrier protein (protein_metabolites.xlsx) stay at 0: mMtoPart floors every ODE
        # species at 1 particle, which would leave 1 per form and a negative base form (protein 0 minus the forms)
        zero = []
        for locus, mode in self.ko.items():
            pid = 'P_' + self.spec['genes'][locus]['num']
            if mode == 'full' and pid in self.spec.get('protein_metabolites', {}):
                zero += self.spec['protein_metabolites'][pid]['forms']
        if zero:
            sim_properties['ko_zero_species'] = zero
            self.log.append(('ODE forms held at 0 (knocked-out carrier)', 'floor 1 particle', ', '.join(zero)))
        if SMC in [s for l in self.ko for s in gene_species(self.spec, l)]:
            sim_properties['smc_min'] = 0          # the knockout removes the max(1, ...) floor on loop extruders

    def _edit_tables(self):
        """Edit the cached kinetic_params.xlsx sheets in place; every ODE and CME build copies from this cache."""
        import processes.Rxns_ODE as ODE
        path = self.head + 'input_data/kinetic_params.xlsx'
        rx = self.spec['reactions']
        edits = []                                   # (reaction, param, new value)
        for e in self.pert.get('reaction_parameters') or []:
            p = next(p for p in rx[e['reaction']]['params'] if p.get('key') == e['parameter'])
            edits.append((rx[e['reaction']], p, p['value'] * e['scale'] if 'scale' in e else e['value']))
        for rid in self.pert.get('disabled_reactions') or []:
            r = rx[rid]
            for p in r['params']:
                if r['kind'] == 'ode_mm' and p.get('key') in ('kcatF', 'kcatR'):
                    edits.append((r, p, 0.0))
                elif r['kind'] != 'ode_mm' and isinstance(p['value'], (int, float)) and p['name'] != 'Radius' and p.get('key'):
                    edits.append((r, p, 0.0))
        for r, p, new in edits:
            sheet, idx = p['src']['sheet'], p['src']['row'] - 2
            ODE._read_excel_cached(path, sheet)                              # make sure the sheet is cached
            key = next(k for k in ODE._EXCEL_CACHE if k[0] == os.path.abspath(path) and k[1] == sheet)
            df = ODE._EXCEL_CACHE[key]
            assert df.at[idx, 'Reaction Name'] == r['id'], (r['id'], sheet, idx)
            old = df.at[idx, 'Value']
            df.at[idx, 'Value'] = float(new)
            self.log.append(('%s.%s' % (r['id'], p.get('key') or p['name']), old, float(new)))

    def _wrap_transcription(self, GIP):
        if not any(k[0] == 'transcription' for k in self.gene_scale):
            return
        orig = GIP.TranscriptionRate
        me = self

        def TranscriptionRate(sim_properties, locusTag, rnasequence):
            return orig(sim_properties, locusTag, rnasequence) * me._rate_scale('transcription', locusTag.split('_')[1])
        GIP.TranscriptionRate = TranscriptionRate
        self.log.append(('GIP_rates.TranscriptionRate', 'computed', 'scaled per gene'))

    def _wrap_sim(self, sim):
        me = self
        rate_orig, diff_orig, dist_orig = sim.rateConst, sim.diffusionConst, sim.distributeNumber

        def rateConst(name, value, order, *a, **kw):
            new = value
            if name in me.rc:
                e = me.rc[name]
                new = value * e['scale'] if 'scale' in e else e['value']
            else:
                for rate, suf in PER_GENE_SUFFIX.items():
                    if len(name) == 4 + len(suf) and name.endswith(suf) and name[:4].isdigit():
                        new = value * me._rate_scale(rate, name[:4])
            if new != value:
                me.log.append(('rate constant ' + name, value, new))
            return rate_orig(name, new, order, *a, **kw)

        def diffusionConst(name, value, *a, **kw):
            s = me.dc.get(name, 1.0)
            if _own_rna_dc(name):
                s *= me.dc.get('rna_diff', 1.0)
            if s != 1.0:
                me.log.append(('diffusion ' + name, value, value * s))
            return diff_orig(name, value * s, *a, **kw)

        def distributeNumber(sp, reg, count):
            name, rname = getattr(sp, 'name', str(sp)), getattr(reg, 'name', str(reg))
            new = int(count)
            if name in me.place:
                target = me.place[name]
                new = 0 if target == 0 else _share(target, rname, me.spec['genes'][me.spec['species'][name]['gene']].get('localization'))
            elif name in me.mrna_mean:
                rng = np.random.RandomState(int(name[2:6]))             # separate stream: the global draws stay as they were
                new = int(rng.poisson(me.mrna_mean[name]))
            if new != int(count):
                me.log.append(('initial %s in %s' % (name, rname), int(count), new))
            me._placed[name] = me._placed.get(name, 0) + new
            return dist_orig(sp, reg, new)

        sim.rateConst, sim.diffusionConst, sim.distributeNumber = rateConst, diffusionConst, distributeNumber

    # ------------------------------------------------------------------------------------------ after particles
    def after_particles(self, sim, sim_properties):
        """Promoters, metabolite pools and medium: set after initializeParticles, before any reaction reads them."""
        genes = self.spec['genes']
        prom = sim_properties['promoters']
        for locus, mode in self.ko.items():
            if mode in ('full', 'expression_only'):
                self.log.append(('promoter ' + locus, prom[locus], 0))
                prom[locus] = 0
        for locus, s in self.kd.items():
            self.log.append(('promoter ' + locus, prom[locus], prom[locus] * s))
            prom[locus] = prom[locus] * s
        from processes.ImportInitialConditions import mMtoPart
        for sid, v in (self.pert.get('initial_metabolites_mM') or {}).items():
            new = mMtoPart(v, sim_properties)
            self.log.append(('initial %s (particles)' % sid, sim_properties['counts'][sid], new))
            sim_properties['counts'][sid] = new
        for sid, v in (self.pert.get('medium_mM') or {}).items():
            self.log.append(('medium %s (mM)' % sid, sim_properties['medium'][sid], v))
            sim_properties['medium'][sid] = v

    # ------------------------------------------------------------------------------------------ read back
    def report(self, sim, sim_properties):
        """Read every edited value back from the simulation objects; write perturbation_applied.json in the run directory."""
        import utility.GIP_rates as GIP
        import processes.Rxns_ODE as ODE
        checks = []

        def check(what, expected, got, tol=1e-9):
            ok = (expected == got) if not isinstance(expected, float) and not isinstance(got, float) else \
                abs(float(expected) - float(got)) <= tol * max(1.0, abs(float(expected)))
            checks.append({'what': what, 'expected': expected, 'got': got, 'ok': bool(ok)})
        # placements: jLM keeps the last count per (region, species)
        dist = dict(self._placed)                                # what the wrapper handed over
        if hasattr(sim, '_particleDistCount'):                  # what jLM holds (last count per region and species)
            names = {s.idx: s.name for s in sim.speciesList}
            dist = {}
            for reg, spi, n in sim._particleDistCount:
                nm = names.get(spi, spi)
                dist[nm] = dist.get(nm, 0) + n
        for sid, target in self.place.items():
            want = 0 if target == 0 else sum(_share(target, r, self.spec['genes'][self.spec['species'][sid]['gene']].get('localization'))
                                             for r in _regions(self.spec['genes'][self.spec['species'][sid]['gene']].get('localization')))
            check('placed %s' % sid, want, dist.get(sid, 0))
        prom = sim_properties.get('promoters', {})
        for locus, mode in self.ko.items():
            if mode in ('full', 'expression_only'):
                check('promoter %s' % locus, 0, prom.get(locus))
        for name, e in self.rc.items():
            obj = sim.rxnRateList[name]
            base = self.spec['rdme']['rate_constants'][name]['value']
            check('rate constant %s' % name, base * e['scale'] if 'scale' in e else e['value'], float(obj.value))
        for (rate, num), s in self.gene_scale.items():
            if rate == 'transcription':
                continue
            nums = [g['num'] for g in self.spec['genes'].values()] if num == '*' else [num]
            for n in nums:
                nm = n + PER_GENE_SUFFIX[rate]
                try:
                    obj = sim.rxnRateList[nm]
                except KeyError:
                    continue
                pg = self.spec['rdme']['per_gene_constants']['<n>' + PER_GENE_SUFFIX[rate]]
                base = next((i['value'] for i in pg['instances'] if i['name'] == nm), None)
                if base is not None:
                    check('rate constant %s' % nm, base * self._rate_scale(rate, n), float(obj.value), 1e-6)
        for name, s in self.dc.items():
            if name == 'rna_diff':
                continue
            check('diffusion %s' % name, self.spec['rdme']['diffusion_constants'][name]['value'] * s, float(sim.diffRateList[name].value))
        for e in self.pert.get('gip_constants') or []:
            check('GIP_rates.%s' % e['constant'], e['value'], getattr(GIP, e['constant']))
        path = self.head + 'input_data/kinetic_params.xlsx'
        for (what, old, new) in [l for l in self.log if '.' in l[0] and not l[0].startswith(('GIP', 'rate', 'diffusion', 'initial', 'medium', 'promoter'))]:
            rid, key = what.split('.', 1)
            r = self.spec['reactions'][rid]
            p = next(p for p in r['params'] if (p.get('key') or p['name']) == key)
            df = ODE._read_excel_cached(path, p['src']['sheet'])
            check('table %s' % what, new, float(df.at[p['src']['row'] - 2, 'Value']))
        for sid, v in (self.pert.get('medium_mM') or {}).items():
            check('medium %s' % sid, v, sim_properties['medium'][sid])
        out = OrderedDict([('name', self.pert.get('name')), ('source', self.path), ('fingerprint', self.spec['meta']['fingerprint']),
                           ('checks_passed', sum(c['ok'] for c in checks)), ('checks_failed', sum(not c['ok'] for c in checks)),
                           ('checks', checks), ('log', [{'what': a, 'from': b, 'to': c} for a, b, c in self.log]),
                           ('resolution', self.resolution['edits'])])
        wd = sim_properties.get('working_directory')
        if wd:
            with open(os.path.join(wd, 'perturbation_applied.json'), 'w') as f:
                json.dump(out, f, indent=1, default=str)
            with open(os.path.join(wd, 'perturbation.yaml'), 'w') as f:
                f.write(P.dump(self.pert))
        print('PERTURBATION %s applied: %d edits logged, %d/%d read-back checks passed' % (
            self.pert.get('name'), len(self.log), out['checks_passed'], len(checks)))
        for c in checks:
            if not c['ok']:
                print('PERTURBATION CHECK FAILED: %s expected %s got %s' % (c['what'], c['expected'], c['got']))
        return out


def _regions(localization):
    if localization == 'trans-membrane':
        return ['membrane']
    if localization == 'peripheral membrane':
        return ['outer_cytoplasm']
    return ['cytoplasm', 'outer_cytoplasm', 'DNA']


def _share(total, region, localization):
    """A requested count split over initializeProteins' regions (1/6 outer cytoplasm, 1/3 DNA, the rest cytoplasm), so exactly
    `total` particles are placed (initializeProteins' own int(n/2)+int(n/6)+int(n/3) can lose one or two to rounding)."""
    total = int(total)
    if localization in ('trans-membrane', 'peripheral membrane'):
        return total
    outer, dna = int(total / 6), int(total / 3)
    return {'cytoplasm': total - outer - dna, 'outer_cytoplasm': outer, 'DNA': dna}.get(region, total)


def load_spec(head):
    """The built spec when it matches the current inputs, otherwise a fresh build (no code index)."""
    from .load import fingerprint, build_spec
    p = os.path.join(head, 'modelspec', 'out', 'model_spec.json')
    fp, _ = fingerprint(head.rstrip('/'))
    if os.path.exists(p):
        spec = json.load(open(p))
        if spec['meta']['fingerprint'] == fp and spec.get('rdme'):
            return spec
        print('PERTURBATION: model_spec.json is stale (inputs changed); rebuilding the spec')
    return build_spec(head.rstrip('/'), with_code=False)
