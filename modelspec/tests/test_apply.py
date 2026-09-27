"""Deterministic gate for modelspec.apply, run in the 4DWCM container from the repo root:

    python -m modelspec.tests.test_apply

Each scenario builds the model on the recording stand-in for the jLM Sim (as modelspec.rdme does) in a fresh process,
seeded, and snapshots everything that enters the model: rate and diffusion constants, t=0 placements, promoters,
metabolite and medium pools, the ODE reaction parameters, the CME rates and the GIP constants.
  1. no perturbation vs an empty perturbation file: identical snapshots
  2. no perturbation vs each test file: exactly the expected entries differ, and nothing else
"""

import json
import os
import subprocess
import sys
import tempfile

HEAD = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, HEAD)

SCENARIOS = {
    'empty': {},
    'combo': {
        'knockouts': [{'gene': 'JCVISYN3A_0227', 'mode': 'full'}, {'gene': 'JCVISYN3A_0415', 'mode': 'full'},
                      {'gene': 'JCVISYN3A_0719', 'mode': 'full'}, {'gene': 'JCVISYN3A_0732', 'mode': 'expression_only'},
                      {'gene': 'JCVISYN3A_0001', 'mode': 'initial_only'}],
        'knockdowns': [{'gene': 'JCVISYN3A_0645', 'promoter_scale': 0.25}],
        'initial_protein_counts': {'P_0609': 5, 'P_0652': 40},
        'initial_mrna_means': {'R_0445': 3.0},
        'initial_metabolites_mM': {'M_atp_c': 1.0},
        'medium_mM': {'M_glc__D_e': 0.0},
        'reaction_parameters': [{'reaction': 'PGI', 'parameter': 'kcatF', 'scale': 0.1},
                                {'reaction': 'PFK', 'parameter': 'Km:M_atp_c', 'value': 0.5},
                                {'reaction': 'ALATRS', 'parameter': 'k_cat', 'value': 15.0},
                                {'reaction': 'GLCpts0', 'parameter': 'kcatF', 'scale': 2}],
        'disabled_reactions': ['LDH_L', 'ACt', 'SERTRS'],
        'rdme_rate_constants': [{'constant': 'ribo_bind', 'scale': 0.5}, {'constant': 'RiboS16Binding', 'value': 1e6}],
        'gene_rate_scales': [{'gene': 'JCVISYN3A_0002', 'rate': 'translation', 'scale': 0.5}, {'gene': '*', 'rate': 'mrna_degradation', 'scale': 2},
                             {'gene': 'JCVISYN3A_0005', 'rate': 'secy_insertion', 'scale': 3}, {'gene': 'JCVISYN3A_0010', 'rate': 'rnap_binding', 'scale': 0.1},
                             {'gene': 'JCVISYN3A_0012', 'rate': 'transcription', 'scale': 4}],
        'diffusion_scales': [{'constant': 'diffPtn', 'scale': 0.5}, {'constant': 'rna_diff', 'scale': 0.1}],
    },
    'gip': {'gip_constants': [{'constant': 'riboKcat', 'value': 6}]},
    'trsc': {'gene_rate_scales': [{'gene': 'JCVISYN3A_0012', 'rate': 'transcription', 'scale': 4}]},
}


# ------------------------------------------------------------------------------------------------ one scenario (child)
def snapshot(pert_path):
    import numpy as np
    np.random.seed(12345)
    from modelspec.rdme import RecordingSim, RecordingCME, REGIONS
    import processes.Diffusion as Diff
    import processes.ImportInitialConditions as IC
    import processes.Rxns_RDME as RDME
    import processes.Rxns_CME as CME
    import processes.Rxns_ODE as ODE
    import processes.MC_RDME_initialization as MC
    import processes.RegionsAndComplexes as RC
    import utility.GIP_rates as GIP
    from Bio import SeqIO
    from modelspec.load import VOLUME_L

    genome = next(SeqIO.parse(HEAD + '/input_data/syn3A.gb', 'gb'))
    dnamap, _ = RC.mapDNA(genome)
    sp = {'head_directory': HEAD + '/', 'genome': dnamap, 'counts': {}, 'rnap_spacing': 400, 'volume_L': VOLUME_L,
          'chromosome_features': {'oriC': {'index': 0}}, 'RiboPtnMap': RC.createRibosomalProteinMap(genome), 'working_directory': None}
    sim = RecordingSim()
    pert = None
    if pert_path:
        from modelspec.apply import Perturbation
        pert = Perturbation.from_file(pert_path, HEAD)
        pert.install(sim, sp)
    region_dict = {r: {'shape': np.zeros((2, 2, 2), bool)} for r in REGIONS}
    # ImportInitialConditions.initializeParticles, without the DNA placement (initializeDnaParticles needs the chromosome)
    Diff.defaultDiffusion(sim, region_dict)
    Diff.generalDiffusionConstants(sim)
    IC.initializeProteins(sim, sp)
    IC.initializeMRNA(sim, sp)
    IC.initializeTRNA(sim, sp)
    IC.initializeRRNA(sim, sp)
    IC.initializeRNAP(sim)
    IC.initializeRibosomeParticles(sim, region_dict)
    IC.initializeDegradosomes(sim, region_dict)
    IC.initializePromoterStrengths(sp)
    IC.initializeMedium(sp)
    IC.initializeMetabolites(sp)
    IC.initializeCostCounters(sp)
    IC.initializeProteinMetabolites(sp)
    IC.initializeTrnaCharging(sp)
    IC.initializeLongRnaTracking(sp)
    if pert:
        pert.after_particles(sim, sp)
    sp['name_to_index'] = {n: i + 1 for i, n in enumerate(sim.species_names)}
    RDME.general_reaction_rates(sim)
    MC.addGeneticInformationReactions(sim, sp)
    RDME.replicationInitiation(sim, sp)
    RDME.addRibosomeBiogenesis(sim, sp)
    RDME.addRNAPassembly(sim, sp)
    csim = RecordingCME()
    for locusTag, d in dnamap.items():
        if d['Type'] == 'rRNA' and len(d['RNAsequence']) > 2 * sp['rnap_spacing']:
            CME.transcriptionLong(csim, sp, locusTag, d['RNAsequence'])
        else:
            CME.transcription(csim, sp, locusTag, d['RNAsequence'])
    CME.tRNAcharging(csim, sp)
    rb = ODE._random_binding_static(HEAD + '/input_data/kinetic_params.xlsx', HEAD + '/input_data/Syn3A_updated.xml')
    nrb = ODE._read_excel_cached(HEAD + '/input_data/kinetic_params.xlsx', 'Non-Random-Binding Reactions')
    import processes.SpatialDnaDynamics as SDD
    sp['counts']['P_0415'] = 0
    snap = {
        'rate_consts': {k: v[0] for k, v in sim.rate_consts.items()},
        'diff_consts': dict(sim.diff_consts),
        'placed': dict(sim.placed),
        'promoters': dict(sp['promoters']),
        'metabolites': {k: v for k, v in sp['counts'].items() if k.startswith('M_')},
        'medium': dict(sp['medium']),
        'ode_rb': {x[0]: [float(x[6]), float(x[7]), [float(k) for k in x[8]], [float(k) for k in x[9]]] for x in rb},
        'ode_nrb': {'%s.%s' % (r['Reaction Name'], r['Parameter Type']): str(r['Value']) for _, r in nrb.iterrows()},
        'cme': {'%s>%s' % ('+'.join(a), '+'.join(b)): r for a, b, r in csim.reactions},
        'gip': {k: float(getattr(GIP, k)) for k in ('riboKcat', 'riboKd', 'ctRNAconc', 'rnaPolKcat', 'rnaPolKd', 'ATPconc')},
        'smc_at_zero': SDD._num_smc(sp),
    }
    return snap


# ------------------------------------------------------------------------------------------------ driver (parent)
def run_child(name, pert):
    path = ''
    if pert is not None:
        from modelspec import perturbation as P
        full = P.empty('test_' + name)
        full.update(pert)
        fd, path = tempfile.mkstemp(suffix='.yaml')
        with os.fdopen(fd, 'w') as f:
            f.write(P.dump(full))
    out = subprocess.run([sys.executable, '-m', 'modelspec.tests.test_apply', '--child', path], cwd=HEAD, capture_output=True, text=True)
    if out.returncode != 0:
        print(out.stdout[-3000:], out.stderr[-3000:])
        raise SystemExit('scenario %s failed to build' % name)
    line = [l for l in out.stdout.splitlines() if l.startswith('SNAPSHOT ')][-1]
    return json.loads(line[len('SNAPSHOT '):])


def diff(a, b):
    out = {}
    for sec in a:
        keys = set(a[sec]) | set(b[sec]) if isinstance(a[sec], dict) else {None}
        if not isinstance(a[sec], dict):
            if a[sec] != b[sec]:
                out[sec] = {'value': (a[sec], b[sec])}
            continue
        d = {k: (a[sec].get(k), b[sec].get(k)) for k in keys if a[sec].get(k) != b[sec].get(k)}
        if d:
            out[sec] = d
    return out


FAIL = []


def check(name, ok, detail=''):
    print('%s  %s%s' % ('ok  ' if ok else 'FAIL', name, ('   ' + detail) if detail and not ok else ''))
    if not ok:
        FAIL.append(name)


def expect_combo(base, spec):
    """The entries the combo file must change (everything else must stay)."""
    g = spec['genes']
    exp = {sec: set() for sec in base}
    exp['promoters'] |= {'JCVISYN3A_0227', 'JCVISYN3A_0415', 'JCVISYN3A_0719', 'JCVISYN3A_0732', 'JCVISYN3A_0645'}
    zero = ['P_0227', 'R_0227', 'P_0415', 'R_0415', 'R_0719', 'P_0001', 'R_0001']      # knocked out or depleted at t=0
    exp['placed'] |= {x for x in zero if base['placed'].get(x, 0) > 0} | {'P_0609', 'P_0652'}
    exp['placed_optional'] = {'R_0445'}                  # new Poisson mean 3: the draw may equal the baseline draw
    exp['metabolites'] |= {'M_atp_c'}
    exp['medium'] |= {'M_glc__D_e'}
    exp['ode_rb'] |= {'PGI', 'PFK', 'LDH_L'}
    exp['ode_nrb'] |= {'GLCpts0.kcatF', 'ACt.P_R'}
    exp['rate_consts'] |= {'ribo_bind', 'RiboS16Binding', '0002_translat', '0005_insertion', '0010RNAP_on'}
    exp['rate_consts'] |= {n + 'RNAP_on' for n in ('0227', '0415', '0719', '0732', '0645')}   # promoter 0 / x0.25 sets their binding rate
    exp['rate_consts'] |= {k for k in base['rate_consts'] if k.endswith('_RNAdeg')}
    exp['diff_consts'] |= {'diffPtn'} | {k for k in base['diff_consts'] if (k.endswith('_diff') or k.endswith('_diffDna')) and k.startswith('R') and not k.startswith('RNAP')}
    # CME: transcription of the knocked-down / scaled genes changes; ATP 1 mM changes every transcription rate (it reads NTP pools);
    # the ALATRS / SERTRS charging reactions change
    exp['cme'] = None
    exp['smc_at_zero'] = {'value'}
    return exp


def main():
    if '--child' in sys.argv:
        path = sys.argv[sys.argv.index('--child') + 1] or None
        print('SNAPSHOT ' + json.dumps(snapshot(path)))
        return
    from modelspec.apply import load_spec
    spec = load_spec(HEAD + '/')
    base = run_child('none', None)
    again = run_child('none', None)
    check('two unperturbed builds are identical (the snapshot is deterministic)', diff(base, again) == {}, str(list(diff(base, again))))
    empty = run_child('empty', SCENARIOS['empty'])
    d = diff(base, empty)
    check('empty perturbation file = no perturbation (%d rate constants, %d placements, %d CME rates compared)' % (
        len(base['rate_consts']), len(base['placed']), len(base['cme'])), d == {}, json.dumps({k: list(v)[:3] for k, v in d.items()}))

    combo = run_child('combo', SCENARIOS['combo'])
    d = diff(base, combo)
    exp = expect_combo(base, spec)
    for sec in base:
        got = set(d.get(sec, {}))
        if exp.get(sec) is None:
            continue
        if sec == 'placed':
            got -= exp['placed_optional']
        check('combo: %-12s changes exactly the expected %d' % (sec, len(exp[sec])), got == exp[sec],
              'unexpected %s; missing %s' % (sorted(got - exp[sec])[:6], sorted(exp[sec] - got)[:6]))
    # value-level checks
    c, b = combo, base
    check('knockout: promoter 0 and nothing placed (pdhC)', c['promoters']['JCVISYN3A_0227'] == 0 and c['placed']['P_0227'] == 0 and c['placed']['R_0227'] == 0)
    check('expression_only: promoter 0, protein still placed (deoC)', c['promoters']['JCVISYN3A_0732'] == 0 and c['placed']['P_0732'] == b['placed']['P_0732'])
    check('initial_only: promoter kept, protein and mRNA zero (dnaA)', c['promoters']['JCVISYN3A_0001'] == b['promoters']['JCVISYN3A_0001'] and c['placed']['P_0001'] == 0)
    check('tRNA knockout places no tRNA (R_0719)', c['placed']['R_0719'] == 0)
    check('knockdown x0.25 (rpoA promoter)', abs(c['promoters']['JCVISYN3A_0645'] - 0.25 * b['promoters']['JCVISYN3A_0645']) < 1e-9)
    check('RNAP binding rate follows the knockdown', abs(c['rate_consts']['0645RNAP_on'] - 0.25 * b['rate_consts']['0645RNAP_on']) < 1e-6 * b['rate_consts']['0645RNAP_on'])
    check('knocked-out gene: RNAP binding rate 0', c['rate_consts']['0227RNAP_on'] == 0)
    check('initial protein count 5 is placed exactly (P_0609)', c['placed']['P_0609'] == 5, 'placed %s' % c['placed']['P_0609'])
    check('initial protein count 40 is placed exactly (P_0652)', c['placed']['P_0652'] == 40, 'placed %s' % c['placed']['P_0652'])
    check('SMC knockout removes the loop-extruder floor', c['smc_at_zero'] == 0 and b['smc_at_zero'] == 1)
    check('PGI kcatF x0.1', abs(c['ode_rb']['PGI'][0] - 0.1 * b['ode_rb']['PGI'][0]) < 1e-9)
    check('PFK Km(ATP) = 0.5', 0.5 in c['ode_rb']['PFK'][2])
    check('LDH_L disabled: kcatF = kcatR = 0', c['ode_rb']['LDH_L'][:2] == [0.0, 0.0])
    check('ACt disabled: P_R = 0', float(c['ode_nrb']['ACt.P_R']) == 0.0)
    check('medium glucose 0', c['medium']['M_glc__D_e'] == 0.0)
    check('ALATRS k_cat = 15 in the CME', any(abs(c['cme'][k] - 15.0) < 1e-9 for k in c['cme'] if k.startswith('P_0163_atp_aa_')))
    check('SERTRS disabled: its CME rates are 0', all(c['cme'][k] == 0 for k in c['cme'] if k.startswith('P_0061')))
    check('CME transcription of the knocked-out gene falls to the floor rate', c['cme']['RP_0227_C1>RP_0227_f_C1'] < b['cme']['RP_0227_C1>RP_0227_f_C1'])
    check('rna_diff x0.1 reaches mRNA and ribosome intermediates', abs(c['diff_consts']['R_0001_diff'] - 0.1 * b['diff_consts']['R_0001_diff']) < 1e-24
          and abs(c['diff_consts']['Rs4s15_diff'] - 0.1 * b['diff_consts']['Rs4s15_diff']) < 1e-24)

    t = run_child('trsc', SCENARIOS['trsc'])
    d = diff(base, t)
    r12 = {k for k in base['cme'] if k.startswith('RP_0012_C')}
    check('transcription x4 for gene 0012: exactly its 2 CME reactions change, each x4', set(d) == {'cme'} and set(d['cme']) == r12 and all(
        abs(t['cme'][k] / base['cme'][k] - 4) < 1e-12 for k in r12), str({k: len(v) for k, v in d.items()}))
    gip = run_child('gip', SCENARIOS['gip'])
    d = diff(base, gip)
    trans = {k for k in base['rate_consts'] if k.endswith('_translat')}
    check('GIP riboKcat=6 changes exactly the %d translation rates (and the GIP value)' % len(trans),
          set(d.get('rate_consts', {})) == trans and set(d) == {'rate_consts', 'gip'}, str(set(d)))
    check('translation rates halve with riboKcat 12 -> 6', all(abs(gip['rate_consts'][k] - 0.5 * base['rate_consts'][k]) < 1e-12 for k in trans))
    print('\n%d failed' % len(FAIL) if FAIL else '\nall checks passed')
    sys.exit(1 if FAIL else 0)


if __name__ == '__main__':
    main()
