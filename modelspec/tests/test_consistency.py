"""ModelSpec against the simulator's own functions. Run in the 4DWCM container from the repo root:

    python -m modelspec.tests.test_consistency

The RDME builders are called on a recording stand-in for the jLM Sim, so what they would add to the lattice (species,
reactions) is compared with the spec without a GPU or a lattice. Exit status 0 = every check passed.
"""

import os
import sys

HEAD = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, HEAD)

FAIL = []


def check(name, ok, detail=''):
    print('%s  %s%s' % ('ok  ' if ok else 'FAIL', name, ('   ' + detail) if detail and not ok else ''))
    if not ok:
        FAIL.append(name)


class _Any:
    """Accepts any attribute access or call; stands in for regions, rate constants and diffusion objects."""
    def __init__(self, name='any'):
        self.name = name

    def __getattr__(self, k):
        return _Any(self.name + '.' + k)

    def __call__(self, *a, **k):
        return _Any(self.name + '()')


class RecSim(_Any):
    def __init__(self):
        super().__init__('sim')
        self.reactions, self.species_names, self.placed = [], set(), {}
        self.rc = _Any('rc')

    def species(self, name):
        self.species_names.add(name)
        return _Any(name)

    def distributeNumber(self, species, region, n):
        self.placed[species.name] = self.placed.get(species.name, 0) + n

    def region(self, name):
        sim = self

        class R(_Any):
            def addReaction(self, subs, prods, rate):
                sim.reactions.append((name, [s.name for s in subs], [p.name for p in prods]))
        return R('region:' + name)


def main():
    from modelspec.load import build_spec
    spec = build_spec(HEAD, with_code=True)
    genes = spec['genes']

    import pandas as pd
    from Bio import SeqIO
    import processes.RegionsAndComplexes as RC
    import processes.ImportInitialConditions as IC
    import processes.Rxns_ODE as ODE
    import processes.Rxns_RDME as RDME

    genome = next(SeqIO.parse(os.path.join(HEAD, 'input_data', 'syn3A.gb'), 'gb'))
    dnamap, _ = RC.mapDNA(genome)
    check('gene set = mapDNA', set(dnamap) == set(genes), str(set(dnamap) ^ set(genes)))
    check('gene types = mapDNA', all(dnamap[k]['Type'] == genes[k]['type'] for k in dnamap))
    sp = {'head_directory': HEAD + '/', 'genome': dnamap, 'counts': {}}

    IC.initializePromoterStrengths(sp)
    bad = [k for k, v in sp['promoters'].items() if abs(v - genes[k]['promoter']) > 1e-9]
    check('promoter strengths = initializePromoterStrengths', not bad, str(bad[:5]))

    sp['volume_L'] = spec['meta']['volume_L']
    IC.initializeMetabolites(sp)
    bad = [k for k, v in sp['counts'].items() if spec['species'][k].get('initial_count') != v]
    check('metabolite particle counts = initializeMetabolites', not bad, str(bad[:5]))
    IC.initializeMedium(sp)
    bad = [k for k, v in sp['medium'].items() if spec['species'][k].get('medium_mM') != v]
    check('medium = initializeMedium', not bad, str(bad[:5]))

    rb = ODE._random_binding_static(HEAD + '/input_data/kinetic_params.xlsx', HEAD + '/input_data/Syn3A_updated.xml')
    bad = []
    for (rid, subs, sst, prods, pst, law, kf, kr, skm, pkm, _rows) in rb:
        r = spec['reactions'][rid]
        if [s for s, _ in r['substrates']] != list(subs) or [-s for _, s in r['substrates']] != list(sst):
            bad.append(rid + ' substrates')
        if [s for s, _ in r['products']] != list(prods) or [s for _, s in r['products']] != list(pst):
            bad.append(rid + ' products')
        num = lambda v: None if v is None else float(v)          # the sheets store many numbers as text
        if (r['kcatF'], r['kcatR']) != (num(kf), num(kr)):
            bad.append(rid + ' kcat')
        kms = {p['species']: p['value'] for p in r['params'] if p.get('key', '').startswith('Km:')}
        if [kms.get(s) for s in subs] != [num(x) for x in skm] or [kms.get(s) for s in prods] != [num(x) for x in pkm]:
            bad.append(rid + ' Km')
    check('%d ODE reactions = _random_binding_static' % len(rb), not bad, str(bad[:5]))
    check('ODE reaction set = spec active ode_mm',
          {x[0] for x in rb} == {k for k, r in spec['reactions'].items() if r['kind'] == 'ode_mm' and r['active']})

    # enzyme rule: getEnzymeConc with one enzyme at a time set to zero must match the spec's requirement groups
    from modelspec.impact import _groups
    sp['volume_L'] = spec['meta']['volume_L']
    bad = []
    for (rid, *_rest, rows) in rb:
        r = spec['reactions'][rid]
        enz = [s for g in r['requires'] if not g.get('via') for s in g['species']]
        for e in enz:
            sp['counts'] = {x: 100 for x in enz}
            sp['counts'][e] = 0
            zero = ODE.getEnzymeConc(rows, sp) == 0
            st, _ = _groups([g for g in r['requires'] if not g.get('via')], {e}, set())
            if zero != (st == 'blocked'):
                bad.append('%s/%s' % (rid, e))
    check('enzyme knockouts = getEnzymeConc', not bad, str(bad[:5]))

    # RDME builders on the recording Sim
    sim = RecSim()
    RDME.addRNAPassembly(sim, sp)
    subs = {s for _, ss, pp in sim.reactions for s in ss if s.startswith('P_')}
    rnap = spec['processes']['rnap_assembly']['requires'][0]['species']
    check('RNAP subunits = addRNAPassembly', subs == set(rnap), str(subs))

    sim = RecSim()
    sp['RiboPtnMap'] = RC.createRibosomalProteinMap(genome)
    RDME.addRibosomeBiogenesis(sim, sp)
    ptns = {s for _, ss, _ in sim.reactions for s in ss if s.startswith('P_')}
    rnas = {s for _, ss, _ in sim.reactions for s in ss if s.startswith('R_') and s[2:].isdigit()}
    ribo = spec['ribosome_assembly']
    check('ribosomal proteins = addRibosomeBiogenesis', ptns == set(ribo['ssu_proteins']) | set(ribo['lsu_proteins']),
          str(ptns ^ (set(ribo['ssu_proteins']) | set(ribo['lsu_proteins']))))
    check('rRNAs in assembly', rnas == {'R_0067', 'R_0068', 'R_0069', 'R_0532', 'R_0533', 'R_0534'}, str(rnas))

    sim = RecSim()
    sp['chromosome_features'] = {'oriC': {'index': 0}}
    RDME.replicationInitiation(sim, sp)
    ptns = {s for _, ss, _ in sim.reactions for s in ss if s.startswith('P_')}
    check('replication proteins = replicationInitiation', ptns == set(spec['processes']['replication_initiation']['requires'][0]['species']), str(ptns))

    sim = RecSim()
    IC.initializeDegradosomes(sim, {'outer_cytoplasm': {'shape': __import__('numpy').zeros((2, 2, 2), bool)}})
    ptns = {s for _, ss, _ in sim.reactions for s in ss if s.startswith('P_')}
    check('degradosome subunits = initializeDegradosomes', ptns == set(spec['processes']['degradosome_assembly']['requires'][0]['species']), str(ptns))

    # membrane proteins: initializeProteins creates C_P_ for exactly the spec's membrane_insertion genes
    sim = RecSim()
    sp['head_directory'] = HEAD + '/'
    IC.initializeProteins(sim, sp)
    cp = {s[2:] for s in sim.species_names if s.startswith('C_P_')}
    mem = {'P_' + genes[l]['num'] for l in spec['processes']['membrane_insertion']['applies_to']}
    check('SecY clients = initializeProteins C_P_ species', cp == mem, str(cp ^ mem))
    bad = [k for k, v in sim.placed.items() if spec['species'][k].get('initial_count') != v]
    check('%d initial protein counts = initializeProteins' % len(sim.placed), not bad and len(sim.placed) == 455, str(bad[:5]))

    # RDME / CME record: every recorded reaction is in exactly one family, and the counts add up
    rd = spec.get('rdme')
    check('RDME recorded in the spec', bool(rd))
    if rd:
        n_fam = sum(f['count'] for f in rd['families'].values())
        check('RDME reactions = sum of family counts (%d)' % rd['n_reactions'], n_fam == rd['n_reactions'], '%d vs %d' % (n_fam, rd['n_reactions']))
        n_cme = sum(f['count'] for f in rd['cme_families'].values())
        check('CME transcription reactions = 2 per gene (+ long rRNA chains): %d' % n_cme, n_cme == rd['n_cme_transcription'])
        n_prot = sum(1 for g in genes.values() if g['type'] == 'protein')
        for fid, expect in [('rdme:translation', n_prot), ('rdme:ribosome_binding', n_prot), ('rdme:mrna_degradation', n_prot),
                            ('rdme:degradosome_binding', n_prot), ('rdme:secy_insertion', len(spec['processes']['membrane_insertion']['applies_to']))]:
            check('%s has one instance per gene (%d)' % (fid, expect), rd['families'][fid]['count'] == expect, str(rd['families'][fid]['count']))
        n_rnap = sum(f['count'] for f in rd['families'].values() if f['family'] == 'rnap_binding')
        check('RNAP binding: 2 copies x %d genes' % len(genes), n_rnap == 2 * len(genes), str(n_rnap))
        check('tRNA charging CME reactions = 5 per tRNA + 2 per synthetase (127)', rd['n_cme_charging'] == 127, str(rd['n_cme_charging']))
        fams = list(rd['families'].values()) + list(rd['cme_families'].values())
        orphans = [f['id'] for f in fams if f.get('process') not in spec['processes']]
        check('every RDME/CME reaction family belongs to a process', not orphans, str(orphans[:5]))
        unassigned = [r for r, x in spec['reactions'].items() if x['active'] and not x.get('process')]
        check('every active spreadsheet reaction belongs to a process', not unassigned, str(unassigned[:5]))
        named = set(rd['rate_constants']) | {i['name'] for pg in rd['per_gene_constants'].values() for i in pg['instances']}
        used = {i['rate_name'] for f in rd['families'].values() for i in f['instances']}
        check('every rate constant an RDME reaction uses is listed (%d)' % len(used), used <= named, str(sorted(used - named)[:5]))
        missing_diff = [row['species'] for f in rd['families'].values() for row in f['diffusion'] if not row.get('constants') and not row.get('immobile')]
        check('every species in an RDME reaction has diffusion or an immobility reason', not missing_diff, str(missing_diff[:5]))
        unused = [k for k, v in rd['rate_constants'].items() if not v['used_by']]
        print('note  rate constants defined but unused: %s' % ', '.join(unused))
        # the recorder's translation rate for a gene equals GIP_rates.TranslationRate on that gene's sequence
        import utility.GIP_rates as GIP
        inst = next(i for i in rd['families']['rdme:translation']['instances'] if i['gene'] == 'JCVISYN3A_0001')
        check('recorded translation rate = GIP.TranslationRate (dnaA)', abs(inst['rate'] - GIP.TranslationRate(dnamap['JCVISYN3A_0001']['AAsequence'])) < 1e-12)
    print('\n%d failed' % len(FAIL) if FAIL else '\nall checks passed')
    sys.exit(1 if FAIL else 0)


if __name__ == '__main__':
    main()
