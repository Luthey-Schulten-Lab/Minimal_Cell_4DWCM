"""Record the RDME and global-CME reactions the simulator builds in Python (nothing in input_data lists them).

The simulator's own builders (ImportInitialConditions, Rxns_RDME, MC_RDME_initialization, Rxns_CME) are run against a
recording stand-in for the jLM Sim, so the list cannot drift from the code. Reactions are grouped into families: the
same reaction written for every gene (translation, degradation, ...) is one family with one instance per gene, and a
one-off reaction (RNAP assembly, replisome loading, ...) is a family of one. Each family carries the rate formula (GIP_rates)
and each instance the value that formula gives for its gene.
"""

import math
import re
from collections import OrderedDict

import numpy as np

REGIONS = ['extracellular', 'membrane', 'outer_cytoplasm', 'cytoplasm', 'ribosomes', 'ribo_centers', 'DNA']


class _Named:
    def __init__(self, name, rec=None):
        self.name, self.rec = name, rec

    def __getattr__(self, k):
        return _Named(self.name + '.' + k, self.rec)

    def __call__(self, *a, **k):
        return self

    def diffusionRate(self, region, dc):                     # species.diffusionRate(region, dc)
        self.rec.transitions.append((self.name, region.name, region.name, dc.name))

    def placeParticle(self, *a):
        pass

    def addReaction(self, subs, prods, rate):                # region.addReaction
        self.rec.reactions.append((self.name, [s.name for s in subs], [p.name for p in prods], rate.name))


class _Table:
    def __init__(self, store, rec):
        self.__dict__['_s'], self.__dict__['_rec'] = store, rec

    def __getattr__(self, k):
        if k not in self._s:
            raise AttributeError('undefined constant %r' % k)
        return _Named(k, self._rec)


class RecordingSim:
    """The subset of jLM.RDME.Sim the builders call, recording instead of building."""

    def __init__(self):
        self.reactions, self.transitions, self.species_names = [], [], []
        self.rate_consts, self.diff_consts = OrderedDict(), OrderedDict()
        self.placed = {}
        self.rc, self.dc = _Table(self.rate_consts, self), _Table(self.diff_consts, self)
        self.diffusionZero = _Named('zero', self)
        self.diff_consts['zero'] = 0.0

    def species(self, name):
        if name not in self.species_names:
            self.species_names.append(name)
        return _Named(name, self)

    def region(self, name):
        return _Named(name, self)

    def rateConst(self, name, value, order):
        self.rate_consts[name] = (float(value), int(order))
        return _Named(name, self)

    def diffusionConst(self, name, value):
        self.diff_consts[name] = float(value)
        return _Named(name, self)

    def transitionRate(self, species, r1, r2, dc):
        self.transitions.append((species.name if species is not None else '*', r1.name, r2.name, dc.name))

    def distributeNumber(self, species, region, n):
        self.placed[species.name] = self.placed.get(species.name, 0) + int(n)


class RecordingCME:
    def __init__(self):
        self.reactions = []

    def addReaction(self, subs, prods, rate):
        tup = lambda x: list(x) if isinstance(x, (list, tuple)) else [x]
        self.reactions.append((tup(subs), tup(prods), float(rate)))


# ------------------------------------------------------------------------------------------------ run the builders
def record(head, genes, ribo_map, volume_L, initial_counts):
    import processes.Diffusion as Diff
    import processes.ImportInitialConditions as IC
    import processes.Rxns_RDME as RDME
    import processes.Rxns_CME as CME
    import processes.MC_RDME_initialization as MC
    import processes.RegionsAndComplexes as RC
    from Bio import SeqIO

    genome = next(SeqIO.parse(head + '/input_data/syn3A.gb', 'gb'))
    dnamap, _ = RC.mapDNA(genome)
    sp = {'head_directory': head + '/', 'genome': dnamap, 'counts': dict(initial_counts), 'rnap_spacing': 400,
          'volume_L': volume_L, 'chromosome_features': {'oriC': {'index': 0}}, 'RiboPtnMap': RC.createRibosomalProteinMap(genome)}
    sim = RecordingSim()
    region_dict = {r: {'shape': np.zeros((2, 2, 2), bool)} for r in REGIONS}
    region_dict['ribo_centers']['shape'] = np.zeros((2, 2, 2), bool)
    np.random.seed(0)                                     # initializeMRNA draws Poisson counts; the draw is not used here
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
    IC.initializeLongRnaTracking(sp)
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
    n_trsc = len(csim.reactions)
    CME.tRNAcharging(csim, sp)
    import utility.GIP_rates as GIP
    consts = {k: getattr(GIP, k) for k in ['rnaPolKcat', 'rnaPolKd', 'rrnaPolKcat', 'riboKcat', 'riboKd', 'ctRNAconc', 'ribosomeConc',
                                            'ATPconc', 'UTPconc', 'CTPconc', 'GTPconc', 'Ecoli_V', 'avgdr', 'countToMiliMol']}
    return sim, csim, n_trsc, consts, sp


# ------------------------------------------------------------------------------------------------ families
NUM = re.compile(r'(?<![0-9])(\d{4})(?![0-9])')            # locus number inside a species or rate name
PER_GENE = {'rnap_binding', 'rnap_binding_next', 'ribosome_binding', 'translation', 'mrna_recycle', 'degradosome_binding',
            'mrna_degradation', 'secy_binding', 'secy_insertion'}
IDX = re.compile(r'_(\d{1,2})(?=_|$)')                  # small index: RNAP position, DnaA site

FAMILY_INFO = {
    # template key -> (id, name, description)
    'decay_removal': ('Decay token removal', 'The "decay" particle left by a degraded mRNA is removed at a very fast rate.'),
}


def _template(name, locus_num):
    """Replace the gene's locus number with <n> and small positional indices with <i>."""
    t = re.sub(r'(?<![0-9])' + locus_num + r'(?![0-9])', '<n>', name) if locus_num else name
    return IDX.sub('_<i>', t)


def _locus_of(names):
    for n in names:
        m = NUM.search(n)
        if m:
            return m.group(1)
    return None


def _describe(subs, prods, region, rate):
    s, p = ' '.join(subs), ' '.join(prods)
    if 'RNAP' in subs and any(x.startswith('G_') for x in subs):
        return 'rnap_binding', 'Transcription: RNAP binds the gene', 'RNA polymerase binds a gene copy on the DNA region; elongation then runs in the global CME.'
    if 'RNAP' in subs and any('_t_' in x for x in subs):
        return 'rnap_binding_next', 'Transcription (rRNA): next RNAP loads behind the previous one', 'Long rRNA genes carry several RNAPs at once (one per 400 nt); each new one binds once the previous has moved on.'
    if 'ribosomeP' in subs and any(x.startswith('R_') for x in subs):
        return 'ribosome_binding', 'Translation: ribosome binds mRNA', 'A free ribosome particle binds an mRNA on a ribosome site; the pair becomes RB_<n>.'
    if any(x.startswith('RB_') for x in subs):
        return 'translation', 'Translation: protein released', 'The ribosome finishes the protein (or its C_P_ precursor for membrane proteins), frees the mRNA as R_<n>_d and drops a cost token.'
    if any(x.endswith('_d') for x in subs):
        return 'mrna_recycle', 'Translation: read mRNA becomes free mRNA again', 'R_<n>_d converts back to R_<n> at a very fast rate, so it can be translated again.'
    if 'Degradosome' in subs and any(x.startswith('R_') for x in subs):
        return 'degradosome_binding', 'mRNA degradation: degradosome binds mRNA', 'A degradosome captures an mRNA in the outer cytoplasm.'
    if any(x.startswith('D_') for x in subs):
        return 'mrna_degradation', 'mRNA degradation: mRNA destroyed, degradosome freed', 'The bound mRNA is degraded at a length-dependent rate; its nucleotides are credited back through the cost accounting.'
    if any(x.startswith('C_P_') for x in subs):
        return 'secy_binding', 'SecY insertion: precursor binds SecY', 'A cytoplasmic membrane-protein precursor C_P_<n> binds SecY (P_0652) in the outer cytoplasm.'
    if any(x.startswith('S_') for x in subs):
        return 'secy_insertion', 'SecY insertion: protein inserted into the membrane', 'SecY releases the protein as the membrane form P_<n>, at a length-dependent rate.'
    if subs == ['decay']:
        return 'decay_removal', 'Decay token removal', 'The token left by a degraded mRNA is removed.'
    if 'P_0001' in subs and 'oriC' in subs:
        return 'dnaa_oric', 'Replication initiation: DnaA binds oriC (high affinity)', 'The first DnaA binds the origin.'
    if 'P_0001' in subs and any(x.startswith('ori_HA') or x.startswith('ori_LA1') for x in subs):
        return 'dnaa_low_affinity', 'Replication initiation: DnaA fills the low-affinity sites', 'Second and third DnaA bind after the high-affinity site is occupied.'
    if 'P_0001' in subs and any(x.startswith('ori_ss') or x.startswith('ori_LA2') for x in subs):
        return 'dnaa_ssdna_on', 'Replication initiation: DnaA binds single-stranded DNA', 'DnaA filaments along 30 ssDNA sites at the unwound origin (site i needs site i-1 filled).'
    if any(x.startswith('ori_ss') for x in subs) and 'P_0001' in prods:
        return 'dnaa_ssdna_off', 'Replication initiation: DnaA unbinds from ssDNA', 'Reverse of the ssDNA binding step.'
    if 'P_0609' in subs and any(x.startswith('ori_rep_') for x in subs):
        return 'replisome_second', 'Replication initiation: second replisome loads', 'The second replisome loads after the first, at sites 20-30.'
    if 'P_0609' in subs:
        return 'replisome_first', 'Replication initiation: first replisome loads', 'Once 20-30 DnaA are bound, a replisome (P_0609) loads on the origin.'
    if prods == ['RNAP'] or any(x.startswith('RNAP_a') for x in prods):
        return 'rnap_assembly', 'RNA polymerase assembly', 'alpha + alpha -> alpha2, + beta -> alpha2-beta, + beta\' -> RNAP.'
    if prods == ['Degradosome']:
        return 'degradosome_assembly', 'Degradosome assembly', 'RNase Y + RNase J1 -> degradosome, in the outer cytoplasm.'
    if prods == ['ribosomeP']:
        return 'ribosome_joining', 'Ribosome formation: 30S + 50S', 'Assembled small and large subunits join into a ribosome.'
    if any(x.startswith('Rs') for x in prods):
        return 'ssu_assembly', '30S assembly step', 'A small-subunit ribosomal protein binds the growing 30S intermediate (reduced Mulder network).'
    if any(x.startswith('R5S') or x.startswith('RL') for x in prods):
        return 'lsu_assembly', '50S assembly step', 'A large-subunit ribosomal protein (or 5S rRNA) binds the growing 50S intermediate (LargeSubunit.xlsx order).'
    return 'other', 'Other RDME reaction', ''


RATE_LAWS = {
    # family id -> (rate-name template, tex, legend rows [(symbol, meaning, link-or-None)], note)
    'translation': ('<n>_translat',
                    r'k_{\mathrm{transl}} = \frac{k_{cat}^{ribo}}{\left(\dfrac{K_d^{ribo}}{[\mathrm{tRNA}]}\right)^{2} + \sum_{a} N_a \dfrac{K_d^{ribo}}{[\mathrm{tRNA}]} + n_{aa} - 1}',
                    [('k_{cat}^{ribo}', 'ribosome elongation rate, riboKcat', 'gip:riboKcat'), ('K_d^{ribo}', 'ribosome-tRNA dissociation constant, riboKd', 'gip:riboKd'),
                     ('[tRNA]', 'charged-tRNA concentration, fixed at ctRNAconc (not the simulated tRNA counts)', 'gip:ctRNAconc'),
                     ('N_a', 'count of amino acid a in the protein', None), ('n_{aa}', 'protein length in amino acids', None)],
                    'GIP_rates.TranslationRate: depends only on the protein sequence. First order in RB_<n>, 1/s.'),
    'rnap_binding': ('<n>RNAP_on',
                     r'k_{on}^{RNAP} = 10 \cdot \frac{180}{765} \cdot \frac{s_{prom}}{180} \cdot \frac{V_{E.coli}\, N_A}{11400 \cdot 60}',
                     [('s_{prom}', 'promoter strength of the gene (proxy from the proteomics count)', None), ('V_{E.coli}', '1e-15 L, Ecoli_V', 'gip:Ecoli_V'),
                      ('N_A', "Avogadro's number", 'gip:avgdr')],
                     'GIP_rates.RNAP_binding: second order (gene + RNAP), 1/M/s. The same value is used for both chromosome copies.'),
    'rnap_binding_next': ('<n>RNAP_on', None, [], 'Same rate constant as the first RNAP binding of that gene.'),
    'mrna_degradation': ('<n>_RNAdeg', r'k_{deg} = \frac{88\ \mathrm{nt/s}}{n_{nt}}',
                         [('n_{nt}', 'transcript length in nucleotides', None)], 'GIP_rates.mrnaDegradationRate: first order in D_<n>, 1/s.'),
    'secy_insertion': ('<n>_insertion', r'k_{ins} = \frac{50}{n_{aa}}',
                       [('n_{aa}', 'protein length in amino acids', None)], 'GIP_rates.TranslocationRate: first order in S_<n>, 1/s.'),
    'cme_transcription': ('transcription',
                          r'k_{trsc} = \frac{\min(85,\ \max(10,\ k_{cat}^{RNAP}\, s_{prom}/180))}{\dfrac{(K_d^{RNAP})^{2}}{C_1 C_2} + \sum_{b} N_b \dfrac{K_d^{RNAP}}{[\mathrm{NTP}_b]} + n_{nt} - 1}',
                          [('k_{cat}^{RNAP}', 'RNAP elongation rate, rnaPolKcat', 'gip:rnaPolKcat'), ('K_d^{RNAP}', 'RNAP-NTP dissociation constant, rnaPolKd', 'gip:rnaPolKd'),
                           ('s_{prom}', 'promoter strength of the gene', None), ('C_1, C_2', 'fixed NTP concentrations for the first two bases (ATPconc, UTPconc, GTPconc, CTPconc)', 'gip:ATPconc'),
                           ('[NTP_b]', 'current concentration of ATP, CTP, GTP or UTP from the ODE (recomputed every hook)', 'M_atp_c'),
                           ('N_b', 'count of base b in the transcript', None), ('n_{nt}', 'transcript length', None)],
                          'GIP_rates.TranscriptionRate, applied in the global CME to RP_<n>_C1/C2 -> RP_<n>_f. The value shown is at the initial NTP concentrations.'),
}

# reaction family -> process (roles.PROCESSES id)
FAMILY_PROCESS = {
    'rnap_binding': 'transcription', 'rnap_binding_next': 'transcription', 'cme_transcription': 'transcription',
    'cme_transcription_long': 'transcription',
    'ribosome_binding': 'translation', 'translation': 'translation', 'mrna_recycle': 'translation',
    'degradosome_binding': 'mrna_degradation', 'mrna_degradation': 'mrna_degradation', 'decay_removal': 'mrna_degradation',
    'secy_binding': 'membrane_insertion', 'secy_insertion': 'membrane_insertion',
    'dnaa_oric': 'replication_initiation', 'dnaa_low_affinity': 'replication_initiation', 'dnaa_ssdna_on': 'replication_initiation',
    'dnaa_ssdna_off': 'replication_initiation', 'replisome_first': 'replication_initiation', 'replisome_second': 'replication_initiation',
    'rnap_assembly': 'rnap_assembly', 'degradosome_assembly': 'degradosome_assembly', 'ssu_assembly': 'ssu_assembly',
    'lsu_assembly': 'lsu_assembly', 'ribosome_joining': 'ribosome_joining',
}

CONSTANT_NOTES = {
    'RNAP_off': 'RNAP unbinding from a gene (not used by any reaction in this build)', 'ribo_bind': '40 * Ecoli_V * N_A / 60 / 6800: ribosome-mRNA association, 1/M/s',
    'Ribo_off': 'ribosome unbinding (not used by any reaction in this build)', 'deg_bind_rate': '11 * N_A * Ecoli_V / 60 / 7800: degradosome-mRNA association, 1/M/s',
    'secY_on': 'SecY-precursor association, 1/M/s', 'secY_off': 'SecY unbinding (not used by any reaction in this build)',
    'RNAP_on': 'generic RNAP binding at promoter strength 180 (each gene uses its own <n>RNAP_on)', 'conversion': 'very fast conversion, 1/s',
    'dnaa_ds_ha': 'DnaA binding the high-affinity oriC site, 1/M/s', 'dnaa_ds_la': 'DnaA binding a low-affinity site, 1/M/s',
    'dnaa_ss_on': 'DnaA binding ssDNA, 1/M/s', 'dnaa_ss_off': 'DnaA leaving ssDNA, 1/s', 'repOn': 'replisome loading, 1/M/s',
    'ptnAssoc': 'RNAP subunit association, 1/M/s', 'degAssoc': 'RNase Y + RNase J1 association, 1/M/s', 'LSU_SSU_bind': '30S + 50S joining, 1/M/s',
}


def build(head, genes, ribo_map, volume_L, initial_counts, species):
    sim, csim, n_trsc, consts, sp = record(head, genes, ribo_map, volume_L, initial_counts)
    fam = OrderedDict()
    seen = set()
    for region, subs, prods, rate in sim.reactions:
        fid, name, desc = _describe(subs, prods, region, rate)
        locus = _locus_of([x for x in subs + prods if x[:2] in ('G_', 'R_', 'RB', 'RP', 'D_', 'DT', 'S_', 'C_')] + [rate]) if fid in PER_GENE else None
        if fid in PER_GENE:
            key = (tuple(_template(s, locus) for s in subs), tuple(_template(p, locus) for p in prods), _template(rate, locus))
        elif fid.startswith('dnaa_ss') or fid.startswith('replisome'):
            key = (tuple(IDX.sub('_<i>', s) for s in subs), tuple(IDX.sub('_<i>', p) for p in prods), rate)
        else:
            key = (tuple(subs), tuple(prods), rate)              # one-off reactions: each is its own family
        f = fam.setdefault((fid, key), {'family': fid, 'name': name, 'description': desc, 'template_subs': list(key[0]),
                                        'template_prods': list(key[1]), 'rate_template': key[2], 'regions': [], 'instances': OrderedDict()})
        if region not in f['regions']:
            f['regions'].append(region)
        inst_key = (tuple(subs), tuple(prods), rate)
        if inst_key in seen:
            continue                                                # same reaction in another region
        seen.add(inst_key)
        value, order = sim.rate_consts[rate]
        vars_ = {}
        if locus:
            vars_['n'] = locus
        m = [IDX.search(x) for x in subs + prods]
        idx = sorted({int(x.group(1)) for x in m if x})
        if idx and '<i>' in ' '.join(key[0] + key[1]):
            vars_['i'] = idx[0]
        inst = {'rate': value, 'order': order, 'vars': vars_, 'rate_name': rate,
                'gene': ('JCVISYN3A_' + locus) if locus and ('JCVISYN3A_' + locus) in genes else None}
        if fid not in PER_GENE:
            inst['subs'], inst['prods'] = subs, prods
        f['instances'][inst_key] = inst

    families = OrderedDict()
    counters = {}
    for (fid, key), f in fam.items():
        counters[fid] = counters.get(fid, 0) + 1
        fkey = fid if counters[fid] == 1 else '%s_%d' % (fid, counters[fid])
        insts = list(f['instances'].values())
        f['instances'] = insts
        f['count'] = len(insts)
        f['id'] = 'rdme:' + fkey
        f['layer'] = 'RDME'
        rl = RATE_LAWS.get(fid)
        if rl:
            f['rate_template'], f['rate_tex'], f['legend'], f['rate_note'] = rl[0], rl[1], rl[2], rl[3]
        rates = [i['rate'] for i in insts]
        f['rate_min'], f['rate_max'], f['rate_median'] = min(rates), max(rates), float(np.median(rates))
        f['order'] = insts[0]['order']
        if f['count'] == 1 and not insts[0]['gene']:
            f['single'] = True
        f['per_gene'] = fid in PER_GENE
        f['process'] = FAMILY_PROCESS.get(fid)
        families[f['id']] = f
    for f in families.values():
        t = ' '.join(f['template_subs'] + f['template_prods'])
        if f['per_gene'] and '_C2' in t:
            f['name'] += ' (chromosome copy 2)'
        elif f['per_gene'] and '_C1' in t:
            f['name'] += ' (chromosome copy 1)'
        if f['family'] in ('rnap_binding',) and '<i>' in t:
            f['name'] = f['name'].replace('binds the gene', 'binds a long rRNA gene')
    # disambiguate assembly families and rRNA-specific ones by naming them after their product
    for f in families.values():
        if f['family'] in ('ssu_assembly', 'lsu_assembly', 'rnap_assembly') and f['count'] == 1:
            i = f['instances'][0]
            f['name'] = '%s: %s + %s -> %s' % (f['name'].split(':')[0], i['subs'][0], i['subs'][1], i['prods'][0][:40])

    # CME
    cme_fams = OrderedDict()
    trsc = {'id': 'cme:transcription', 'family': 'cme_transcription', 'layer': 'CME', 'name': 'Transcription elongation (global CME)',
            'description': 'RP_<n>_C1/C2 (RNAP on the gene) -> RP_<n>_f (transcript finished). One reaction per gene copy; the rate follows the '
                           'current NTP pools and is rebuilt every hook.', 'regions': ['global CME'], 'instances': [], 'template_subs': ['RP_<n>_C1'],
            'template_prods': ['RP_<n>_f_C1']}
    rl = RATE_LAWS['cme_transcription']
    trsc['rate_template'], trsc['rate_tex'], trsc['legend'], trsc['rate_note'] = rl
    long_ = {'id': 'cme:transcription_long', 'family': 'cme_transcription_long', 'layer': 'CME', 'name': 'Transcription elongation, multi-RNAP rRNA (global CME)',
             'description': 'For the long rRNA genes each RNAP advances one 400-nt segment at a time (RP_<n>_c<chromosome>_<i> states), so several RNAPs '
                            'transcribe the same gene at once.', 'regions': ['global CME'], 'instances': [], 'rate_template': 'transcription (per segment)',
             'rate_tex': rl[1], 'legend': rl[2], 'rate_note': rl[3], 'template_subs': ['RP_<n>_c1_open_<i+1>', 'RP_<n>_c1_<i>'], 'template_prods': ['RP_<n>_c1_<i+1>', 'RP_<n>_c1_open_<i>']}
    for subs, prods, rate in csim.reactions[:n_trsc]:
        locus = _locus_of(subs + prods)
        inst = {'rate': rate, 'order': 1, 'vars': {'n': locus}, 'gene': 'JCVISYN3A_' + locus, 'rate_name': 'transcription'}
        if len(subs) == 2:
            inst['subs'], inst['prods'] = subs, prods
        (long_ if len(subs) == 2 else trsc)['instances'].append(inst)
    for f in (trsc, long_):
        f['per_gene'] = True
        f['process'] = FAMILY_PROCESS[f['family']]
        rates = [i['rate'] for i in f['instances']]
        f['count'] = len(rates)
        f['rate_min'], f['rate_max'], f['rate_median'] = min(rates), max(rates), float(np.median(rates))
        f['order'] = 1
        cme_fams[f['id']] = f
    n_charging = len(csim.reactions) - n_trsc

    # diffusion: profiles = the set of (from, to, dc) a species has; species -> profile
    # coefficients made by Diffusion.rnaDiffusion (one per RNA-containing species, from its length): '<species>_diff' / '_diffDna'
    def own_dc(dc):
        base = dc[:-8] if dc.endswith('_diffDna') else dc[:-5] if dc.endswith('_diff') else None
        return base is not None and base.startswith('R') and base != 'RNAP'
    prof_of, profiles = {}, OrderedDict()
    per_species = OrderedDict()
    for s, r1, r2, dc in sim.transitions:
        per_species.setdefault(s, []).append((r1, r2, dc))
    for s, lst in per_species.items():
        gen = tuple(sorted({(r1, r2, ('rna_diffDna' if dc.endswith('Dna') else 'rna_diff') if own_dc(dc) else dc) for r1, r2, dc in lst}))
        pid = profiles.get(gen)
        if pid is None:
            pid = 'profile_%d' % (len(profiles) + 1)
            profiles[gen] = pid
        own = {dc: sim.diff_consts[dc] for _, _, dc in lst if own_dc(dc)}
        prof_of[s] = {'profile': pid, 'own': own}
    profiles_by_id = {pid: gen for gen, pid in profiles.items()}
    profile_list = OrderedDict()
    for gen, pid in profiles.items():
        members = [s for s, v in prof_of.items() if v['profile'] == pid]
        kinds = sorted({(species.get(m, {}).get('kind') or ('membrane precursor' if m.startswith('C_P_') else 'other')) for m in members})
        profile_list[pid] = {'id': pid, 'transitions': [list(t) for t in gen], 'members': len(members), 'kinds': kinds,
                             'example': members[0], 'members_list': members if len(members) <= 12 else None}

    def immobile_reason(x):
        if x.startswith(('G_', 'RP_')) or x.startswith('ori'):
            return 'DNA-bound: no diffusion rule, so it only moves when the chromosome model (btree_chromo) repositions the DNA'
        if x.endswith('_TC'):
            return 'translation-cost token: no diffusion rule; read and removed by the hook (Communicate.calculateTranslationCosts)'
        if x.startswith('S_'):
            return 'bound to SecY in the membrane: no diffusion rule, so it stays put until inserted'
        return 'no diffusion rule: the default (zero between all regions) applies, so it does not move'
    fam_diff = {}
    for f in list(families.values()):
        rows = OrderedDict()
        for t in f['template_subs'] + f['template_prods']:
            if t in rows:
                continue
            names = []
            for i in f['instances']:
                if 'subs' in i:
                    names += [x for x in i['subs'] + i['prods'] if _template(x, i['vars'].get('n')) == t or x == t]
                else:
                    x = t.replace('<n>', i['vars'].get('n') or '<n>')
                    if i['vars'].get('i') is not None:
                        x = x.replace('<i+1>', str(i['vars']['i'] + 1)).replace('<i>', str(i['vars']['i']))
                    names.append(x)
            names = list(OrderedDict.fromkeys(names)) or [t]
            known = [n for n in names if n in prof_of]
            if not known:
                rows[t] = {'species': t, 'example': names[0], 'immobile': immobile_reason(names[0])}
                continue
            dcs = OrderedDict()
            for n in known:
                d = prof_of[n]
                for (_, _, dc) in profiles_by_id[d['profile']]:
                    real = dc
                    if dc.startswith('rna_diff'):
                        real = next((k for k in d['own'] if k.endswith('Dna') == dc.endswith('Dna')), dc)
                        val = d['own'].get(real)
                        key = 'RNA length-dependent' + (' (in DNA)' if dc.endswith('Dna') else '')
                    else:
                        val, key = sim.diff_consts.get(dc), dc
                    if key == 'zero' or val is None:
                        continue
                    dcs.setdefault(key, []).append(val)
            rows[t] = {'species': t, 'example': known[0], 'profiles': sorted({prof_of[n]['profile'] for n in known}),
                       'constants': [{'name': k, 'min': min(v), 'max': max(v)} for k, v in dcs.items()]}
        fam_diff[f['id']] = list(rows.values())
    for f in families.values():
        f['diffusion'] = fam_diff.get(f['id'], [])
    immobile = OrderedDict()
    for f in families.values():
        for i in f['instances']:
            for x in i.get('subs', []) + i.get('prods', []):
                if x not in prof_of and x not in immobile:
                    immobile[x] = immobile_reason(x)
    consts_out = OrderedDict((k, {'value': float(v)}) for k, v in consts.items())
    consts_out['ctRNAconc']['note'] = '200 charged tRNA per cell, as mM; TranslationRate uses this constant, not the simulated tRNA counts'
    consts_out['ribosomeConc']['note'] = '500 ribosomes per cell, as mM (not used by the rates)'
    for k, unit in [('rnaPolKcat', 'nt/s'), ('rrnaPolKcat', 'nt/s'), ('riboKcat', 'aa/s'), ('rnaPolKd', 'mM'), ('riboKd', 'mM'), ('ATPconc', 'mM'),
                    ('UTPconc', 'mM'), ('CTPconc', 'mM'), ('GTPconc', 'mM'), ('Ecoli_V', 'L'), ('avgdr', '1/mol'), ('countToMiliMol', 'mM per particle')]:
        consts_out[k]['unit'] = unit
    consts_out['rrnaPolKcat']['note'] = 'defined but the transcription rate caps at 85 nt/s via min(85, ...)'

    used_by = {}
    for f in families.values():
        for i in f['instances']:
            used_by.setdefault(i['rate_name'], [])
            if f['id'] not in used_by[i['rate_name']]:
                used_by[i['rate_name']].append(f['id'])
    rate_consts, per_gene = OrderedDict(), OrderedDict()
    PER_GENE_RC = {'<n>RNAP_on': ('RNAP binding to the gene (both copies)', 'GIP_rates.RNAP_binding: from the promoter strength', 'rnap_binding'),
                   '<n>_translat': ('translation', 'GIP_rates.TranslationRate: from the protein sequence', 'translation'),
                   '<n>_RNAdeg': ('mRNA degradation once bound', 'GIP_rates.mrnaDegradationRate: 88 nt/s / length', 'mrna_degradation'),
                   '<n>_insertion': ('SecY insertion', 'GIP_rates.TranslocationRate: 50 / length', 'secy_insertion')}
    for name, (value, order) in sim.rate_consts.items():
        m = NUM.match(name)
        if m and name[4:] in ('RNAP_on', '_translat', '_RNAdeg', '_insertion'):
            t = '<n>' + name[4:]
            pg = per_gene.setdefault(t, {'template': t, 'what': PER_GENE_RC[t][0], 'formula': PER_GENE_RC[t][1], 'gene_rate': PER_GENE_RC[t][2],
                                         'order': order, 'unit': '1/s' if order == 1 else '1/M/s', 'used_by': [], 'instances': []})
            pg['instances'].append({'name': name, 'gene': 'JCVISYN3A_' + m.group(1), 'value': value})
            for fid in used_by.get(name, []):
                if fid not in pg['used_by']:
                    pg['used_by'].append(fid)
            continue
        group = 'assembly' if (name.startswith('Ribo') and name.endswith('Binding')) or name.endswith('binding') or name == 'LSU_SSU_bind' else 'shared'
        note = CONSTANT_NOTES.get(name, '')
        if not note and name.startswith('Ribo') and name.endswith('Binding'):
            note = '30S assembly: %s binding rate from oneParamMulder-local_min.json (x 1e6), 1/M/s' % name[4:-7]
        elif not note and name.endswith('binding'):
            note = '50S assembly: %s binding rate from LargeSubunit.xlsx (x 1e6), 1/M/s' % name[:-7]
        rate_consts[name] = {'value': value, 'order': order, 'unit': '1/s' if order == 1 else '1/M/s', 'used_by': used_by.get(name, []),
                             'note': note, 'group': group}
    for pg in per_gene.values():
        v = [i['value'] for i in pg['instances']]
        pg['count'], pg['min'], pg['median'], pg['max'] = len(v), min(v), float(np.median(v)), max(v)
    return {'regions': REGIONS, 'families': families, 'cme_families': cme_fams, 'n_reactions': len(seen), 'n_reaction_entries': len(sim.reactions),
            'n_cme_transcription': n_trsc, 'n_cme_charging': n_charging, 'rate_constants': rate_consts, 'per_gene_constants': per_gene,
            'immobile_note': 'Species with no diffusion rule of their own get the default set by Diffusion.defaultDiffusion: zero between all regions.',
            'immobile_patterns': {'G_<n>_C1/C2, RP_<n>_*, oriC, ori_*_DnaA*': immobile_reason('G_'), 'P_<n>_TC': immobile_reason('x_TC'), 'S_<n>': immobile_reason('S_')},
            'diffusion_constants': OrderedDict((k, {'value': v, 'unit': 'm^2/s'}) for k, v in sim.diff_consts.items() if not own_dc(k)),
            'rna_diff_range': [min(v for k, v in sim.diff_consts.items() if own_dc(k) and not k.endswith('Dna')),
                               max(v for k, v in sim.diff_consts.items() if own_dc(k) and not k.endswith('Dna'))],
            'rna_diff_count': sum(1 for k in sim.diff_consts if own_dc(k) and not k.endswith('Dna')),
            'diffusion_profiles': profile_list, 'species_diffusion': prof_of, 'gip_constants': consts_out}
