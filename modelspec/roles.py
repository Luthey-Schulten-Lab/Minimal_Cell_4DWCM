"""Roles the simulator hardcodes in Python rather than reading from input_data/.

Each process names the functions that implement it as 'path::function'; coderefs resolves them to line numbers at build
time and fails the build if a function has been renamed or removed, so this table cannot silently drift from the code.

`requires` is a list of groups: {'rule': 'all'|'any', 'species': [...]}. An 'all' group is broken by losing any member,
an 'any' group only by losing every member (partial loss = reduced). `on_loss` downgrades a broken group ('reduced').
`depends_on` propagates: a process whose dependency is blocked is itself blocked (or reduced, per `effect`).
"""

RNAP_SUBUNITS = ['P_0645', 'P_0804', 'P_0803']            # Rxns_RDME.addRNAPassembly: alpha, beta, beta'
DEGRADOSOME = ['P_0359', 'P_0600']                        # ImportInitialConditions.initializeDegradosomes: RNase Y, RNase J1
DNAA, REPLISOME, SECY, SMC = 'P_0001', 'P_0609', 'P_0652', 'P_0415'
SSU = 'Rs3s4s5s6s7s8s9s10s11s12s13s14s15s16s17s19s20'
LSU = 'R5SL1L2L3L4L5L6L7L9L10L11L13L14L15L16L17L18L19L20L21L22L23L24L27L28L29L31L32L33L34L35L36'

# metabolites the hook reads directly (costs, membrane area); not 'unused' even when no ODE reaction touches them
HOOK_METABOLITES = {'M_atp_c', 'M_ctp_c', 'M_gtp_c', 'M_utp_c', 'M_datp_c', 'M_dttp_c', 'M_dctp_c', 'M_dgtp_c', 'M_amp_c',
                    'M_ump_c', 'M_cmp_c', 'M_gmp_c', 'M_clpn_c', 'M_chsterol_c', 'M_sm_c', 'M_pc_c', 'M_pg_c', 'M_galfur12dgr_c',
                    'M_12dgr_c', 'M_pa_c', 'M_cdpdag_c', 'M_10fthfglu3_c', 'M_thfglu3_c'}

COMPLEXES = {
    'RNAP': {'name': 'RNA polymerase (lumped)', 'initial_count': 0,
             'note': 'no RNAP is placed at t=0 (the distributeNumber in initializeRNAP is commented out): every polymerase is '
                     'assembled from P_0645 x2 + P_0804 + P_0803'},
    'RNAP_aa': {'name': 'RNAP assembly intermediate (alpha2)', 'initial_count': 0},
    'RNAP_aab1': {'name': 'RNAP assembly intermediate (alpha2 beta)', 'initial_count': 0},
    'Degradosome': {'name': 'Degradosome (RNase Y + RNase J1)', 'initial_count': 0,
                    'note': 'none placed at t=0; assembled from P_0359 + P_0600 in the outer cytoplasm'},
    'ribosomeP': {'name': 'Ribosome (70S particle)', 'initial_count': 'from geometry',
                  'note': 'placed at t=0 on every ribo_centers site (initializeRibosomeParticles); new ones come only from '
                          'subunit assembly'},
    SSU: {'name': '30S small subunit (assembled)', 'initial_count': 0},
    LSU: {'name': '50S large subunit (assembled)', 'initial_count': 0},
}

NAMING = [
    ['G_<n>_C1 / _C2', 'gene copy on chromosome 1 / 2 (RDME DNA particle)'],
    ['RP_<n>_C1', 'RNAP bound to the gene (transcribing)'],
    ['R_<n>', 'RNA: mRNA, tRNA or rRNA of locus n'],
    ['R_<n>_d', 'mRNA released by a ribosome (converts back to R_<n>)'],
    ['R_<n>_ch', 'charged tRNA'],
    ['RB_<n>', 'mRNA bound by a ribosome (with _cp, _p, _pe polysome states)'],
    ['P_<n>', 'protein of locus n (for membrane proteins: the inserted form)'],
    ['C_P_<n>', 'cytoplasmic precursor of a trans-membrane protein, before SecY insertion'],
    ['S_<n>', 'membrane protein bound to SecY'],
    ['P_<n>_TC', 'translation-cost token'],
    ['D_<n> / DT_<n>', 'mRNA bound to a degradosome / degradation tracker'],
    ['PM_<n>, RPM_<n>, DM_<n>', 'counters: proteins made, transcripts made, mRNAs degraded'],
    ['M_<id>_c / _e', 'metabolite, cytoplasm / medium'],
]

G, P = 'per_gene', 'global'

PROCESSES = [
    {'id': 'transcription', 'name': 'Transcription', 'layer': 'RDME + CME', 'scope': G,
     'description': 'RNAP binds G_<n> at a rate set by the promoter strength, elongates in the global CME, and releases '
                    'R_<n>. rRNA longer than 800 nt is transcribed by several RNAP at once (transcriptionLong).',
     'anchors': ['processes/Rxns_RDME.py::transcription', 'processes/Rxns_RDME.py::transcriptionLong',
                 'processes/Rxns_CME.py::transcription', 'processes/Rxns_CME.py::transcriptionLong',
                 'utility/GIP_rates.py::RNAP_binding', 'utility/GIP_rates.py::TranscriptionRate',
                 'processes/ImportInitialConditions.py::initializePromoterStrengths',
                 'processes/Communicate.py::updateTranscriptionStates', 'processes/Communicate.py::calculateTranscriptionCosts'],
     'requires': [{'rule': 'all', 'species': ['RNAP']}],
     'depends_on': [{'process': 'rnap_assembly', 'effect': 'blocked',
                     'verified': 'initializeRNAP places no RNAP; addRNAPassembly is the only producer'}],
     'species_roles': {'RNAP': 'polymerase', 'M_atp_c': 'NTP cost', 'M_ctp_c': 'NTP cost', 'M_gtp_c': 'NTP cost',
                       'M_utp_c': 'NTP cost'}},
    {'id': 'translation', 'name': 'Translation', 'layer': 'RDME', 'scope': G,
     'description': 'A ribosome particle binds R_<n> (RB_<n>), then releases the protein (or its C_P_ precursor) and a '
                    'cost token paid in charged tRNA and GTP by the CME and ODE.',
     'anchors': ['processes/Rxns_RDME.py::translation', 'utility/GIP_rates.py::TranslationRate',
                 'processes/MC_RDME_initialization.py::mapTranslationStates', 'processes/RibosomesRDME.py::placeRibosomes',
                 'processes/Communicate.py::calculateTranslationCosts'],
     'requires': [{'rule': 'all', 'species': ['ribosomeP']}],
     'depends_on': [{'process': 'ribosome_joining', 'effect': 'reduced',
                     'verified': 'ribosomeP is only produced by SSU + LSU binding (addRibosomeBiogenesis)'}],
     'species_roles': {'ribosomeP': 'ribosome'},
     'consequence': 'the translation rate uses a fixed tRNA concentration (GIP_rates.TranslationRate), so a blocked tRNA '
                    'synthetase does not slow translation; only the amino-acid cost for that residue is never paid'},
    {'id': 'mrna_degradation', 'name': 'mRNA degradation', 'layer': 'RDME', 'scope': G,
     'description': 'The degradosome binds R_<n> in the outer cytoplasm and returns NMPs to metabolism.',
     'anchors': ['processes/Rxns_RDME.py::degradation_mrna', 'utility/GIP_rates.py::mrnaDegradationRate',
                 'processes/Communicate.py::calculateDegradationCosts'],
     'requires': [{'rule': 'all', 'species': ['Degradosome']}],
     'depends_on': [{'process': 'degradosome_assembly', 'effect': 'blocked',
                     'verified': 'initializeDegradosomes places none; P_0359 + P_0600 is the only producer'}],
     'species_roles': {'Degradosome': 'nuclease complex'}},
    {'id': 'membrane_insertion', 'name': 'SecY membrane insertion', 'layer': 'RDME', 'scope': G,
     'description': 'Trans-membrane proteins are made as C_P_<n> and inserted into the membrane by SecY (P_0652).',
     'anchors': ['processes/Rxns_RDME.py::translocation_secy', 'utility/GIP_rates.py::TranslocationRate'],
     'requires': [{'rule': 'all', 'species': [SECY]}], 'species_roles': {SECY: 'translocon'}},

    {'id': 'rnap_assembly', 'name': 'RNA polymerase assembly', 'layer': 'RDME', 'scope': P,
     'description': 'alpha + alpha -> alpha2; + beta -> alpha2 beta; + beta\' -> RNAP. RNAP starts at zero, so this is the '
                    'only source of polymerase.',
     'anchors': ['processes/Rxns_RDME.py::addRNAPassembly', 'processes/ImportInitialConditions.py::initializeRNAP'],
     'requires': [{'rule': 'all', 'species': RNAP_SUBUNITS}],
     'species_roles': {'P_0645': 'RNAP subunit (alpha, x2)', 'P_0804': 'RNAP subunit (beta)', 'P_0803': "RNAP subunit (beta')",
                       'RNAP': 'product', 'RNAP_aa': 'intermediate', 'RNAP_aab1': 'intermediate'},
     'consequence': 'no polymerase is ever made, so no gene is transcribed'},
    {'id': 'degradosome_assembly', 'name': 'Degradosome assembly', 'layer': 'RDME', 'scope': P,
     'description': 'RNase Y + RNase J1 -> Degradosome in the outer cytoplasm. None exist at t=0.',
     'anchors': ['processes/ImportInitialConditions.py::initializeDegradosomes'],
     'requires': [{'rule': 'all', 'species': DEGRADOSOME}],
     'species_roles': {'P_0359': 'degradosome subunit (RNase Y)', 'P_0600': 'degradosome subunit (RNase J1)',
                       'Degradosome': 'product'},
     'consequence': 'mRNA is never degraded and NMPs are not recycled'},
    {'id': 'replication_initiation', 'name': 'Replication initiation', 'layer': 'RDME', 'scope': P,
     'description': 'DnaA (P_0001) fills oriC high- and low-affinity sites and 30 ssDNA sites; the replisome (P_0609) loads '
                    'at 20-30 bound DnaA. Replication starts once ori_rep2_DnaA_* is occupied.',
     'anchors': ['processes/Rxns_RDME.py::replicationInitiation', 'processes/Communicate.py::checkRepInitState'],
     'requires': [{'rule': 'all', 'species': [DNAA, REPLISOME]}],
     'species_roles': {DNAA: 'initiator (DnaA)', REPLISOME: 'replisome loader'},
     'consequence': 'the chromosome is never replicated; division is still triggered by membrane growth alone'},
    {'id': 'dna_replication', 'name': 'DNA replication (elongation)', 'layer': 'BD + hook', 'scope': P,
     'description': 'After initiation, each DNA hook advances both forks at a dNTP-limited rate and btree_chromo grows the '
                    'daughter chromosome.',
     'anchors': ['processes/SpatialDnaDynamics.py::getReplicatedSegments', 'utility/GIP_rates.py::ReplicationRate',
                 'processes/SpatialDnaDynamics.py::updateChromosome'],
     'requires': [], 'depends_on': [{'process': 'replication_initiation', 'effect': 'blocked',
                                     'verified': 'Hook sets rep_started only when checkRepInitState sees ori_rep2_DnaA_* occupied'}],
     'species_roles': {'M_datp_c': 'dNTP cost', 'M_dttp_c': 'dNTP cost', 'M_dctp_c': 'dNTP cost', 'M_dgtp_c': 'dNTP cost'}},
    {'id': 'ssu_assembly', 'name': '30S subunit assembly', 'layer': 'RDME', 'scope': P,
     'description': '16S rRNA (R_0069 or R_0534) binds small-subunit proteins in the reduced Mulder network '
                    '(oneParamMulder-local_min.json, <= 19 intermediates).',
     'anchors': ['processes/Rxns_RDME.py::addRibosomeBiogenesis'], 'requires': 'ssu', 'species_roles': {SSU: 'product'}},
    {'id': 'lsu_assembly', 'name': '50S subunit assembly', 'layer': 'RDME', 'scope': P,
     'description': '23S rRNA (R_0068 or R_0533) binds 5S rRNA (R_0067 or R_0532) and large-subunit proteins in the order '
                    'of LargeSubunit.xlsx.',
     'anchors': ['processes/Rxns_RDME.py::addRibosomeBiogenesis'], 'requires': 'lsu', 'species_roles': {LSU: 'product'}},
    {'id': 'ribosome_joining', 'name': 'Ribosome formation (30S + 50S)', 'layer': 'RDME', 'scope': P,
     'description': 'SSU + LSU -> ribosomeP. The ribosomes placed at t=0 keep working; this only adds new ones.',
     'anchors': ['processes/Rxns_RDME.py::addRibosomeBiogenesis', 'processes/ImportInitialConditions.py::initializeRibosomeParticles'],
     'requires': [], 'depends_on': [{'process': 'ssu_assembly', 'effect': 'blocked', 'verified': 'SSU + LSU -> ribosomeP is the only joining reaction'},
                                    {'process': 'lsu_assembly', 'effect': 'blocked', 'verified': 'SSU + LSU -> ribosomeP is the only joining reaction'}],
     'species_roles': {'ribosomeP': 'product'},
     'consequence': 'no new ribosomes; translation slows as the initial ribosomes are diluted by growth and division'},
    {'id': 'loop_extrusion', 'name': 'SMC loop extrusion', 'layer': 'BD', 'scope': P,
     'description': 'numSmc for btree_chromo = max(1, int(P_0415 / 2 * dna_smc_bound_fraction)), rewritten every DNA hook.',
     'anchors': ['processes/SpatialDnaDynamics.py::_num_smc', 'processes/SpatialDnaDynamics.py::_write_loop_params_file'],
     'requires': [{'rule': 'all', 'species': [SMC], 'on_loss': 'reduced'}], 'species_roles': {SMC: 'SMC (loop extruder)'},
     'consequence': 'max(1, ...) keeps one looper even at zero SMC protein, so a knockout is not a full loss of looping'},
    {'id': 'trna_charging', 'name': 'tRNA charging', 'layer': 'CME', 'scope': P,
     'description': 'Each synthetase charges the tRNAs of its amino acid (reactions *TRS, tRNA Charging sheet).',
     'anchors': ['processes/Rxns_CME.py::tRNAcharging', 'processes/ImportInitialConditions.py::initializeTrnaCharging',
                 'processes/ImportInitialConditions.py::setCmeSpeciesList'], 'requires': []},
    {'id': 'cost_accounting', 'name': 'Energy and monomer costs', 'layer': 'hook', 'scope': P,
     'description': 'Transcription, translation, degradation and replication accumulate NTP, dNTP and amino-acid costs that '
                    'are paid out of the ODE metabolite pools each hook.',
     'anchors': ['processes/ImportInitialConditions.py::initializeCostCounters',
                 'processes/Communicate.py::communicateCostsToMetabolism', 'processes/Communicate.py::resetCostCounters'],
     'requires': [],
     'species_roles': {'M_atp_c': 'transcription (1 per nt), mRNA degradation, replication, SecY translocation, and ATP in transcripts',
                       'M_gtp_c': 'translation (2 per amino acid) and GTP in transcripts', 'M_ctp_c': 'CTP in transcripts', 'M_utp_c': 'UTP in transcripts',
                       'M_datp_c': 'dATP in replicated DNA', 'M_dttp_c': 'dTTP in replicated DNA', 'M_dctp_c': 'dCTP in replicated DNA',
                       'M_dgtp_c': 'dGTP in replicated DNA', 'M_amp_c': 'returned by mRNA degradation', 'M_ump_c': 'returned by mRNA degradation',
                       'M_cmp_c': 'returned by mRNA degradation', 'M_gmp_c': 'returned by mRNA degradation',
                       'M_10fthfglu3_c': 'formyl donor for fMet (translation initiation cost)', 'M_thfglu3_c': 'returned after formylation'}},
    {'id': 'membrane_growth', 'name': 'Membrane growth', 'layer': 'hook', 'scope': P,
     'description': 'Surface area is summed from lipid counts x headgroup area plus membrane proteins x 28 nm2; the cell '
                    'radius follows it.',
     'anchors': ['processes/Communicate.py::updateSA', 'processes/Growth.py::grow_cell'], 'requires': [],
     'species_roles': {'M_clpn_c': 'cardiolipin, 0.40 nm2 each', 'M_chsterol_c': 'cholesterol, 0.35 nm2', 'M_sm_c': 'sphingomyelin, 0.45 nm2',
                       'M_pc_c': 'phosphatidylcholine, 0.55 nm2', 'M_pg_c': 'phosphatidylglycerol, 0.60 nm2',
                       'M_galfur12dgr_c': 'galactolipid, 0.60 nm2', 'M_12dgr_c': 'diacylglycerol, 0.50 nm2',
                       'M_pa_c': 'phosphatidic acid, 0.50 nm2', 'M_cdpdag_c': 'CDP-diacylglycerol, 0.50 nm2'}},
    {'id': 'division', 'name': 'Cell division', 'layer': 'hook', 'scope': P,
     'description': 'Starts when the surface area implies more than double the initial volume (updateSA sets '
                    'division_started); from then on each DNA hook calls divide_cell instead of grow_cell. The trigger '
                    'is membrane area only: it does not wait for replication or chromosome partitioning.',
     'anchors': ['processes/Communicate.py::updateSA', 'processes/Division.py::divide_cell',
                 'processes/SpatialDnaDynamics.py::partitionChromosomes'],
     'requires': [],
     'consequence': 'a knockout that stops replication does not stop division in the model; the division just '
                    'goes ahead without a replicated chromosome'},
    {'id': 'ode_enzymes', 'name': 'Enzyme concentrations in the ODE', 'layer': 'ODE', 'scope': P,
     'description': 'Each metabolic rate is scaled by its enzyme concentration from the protein count: one enzyme, the sum '
                    'for "or", the minimum for "and", or a fixed 0.001 mM for "default".',
     'anchors': ['processes/Rxns_ODE.py::getEnzymeConc', 'processes/Rxns_ODE.py::defineRandomBindingRxns',
                 'processes/Rxns_ODE.py::defineNonRandomBindingRxns', 'processes/Rxns_ODE.py::addProteinMetabolites'],
     'requires': []},
]


def build_processes(genes, ribo, trna_by_aa):
    out = {}
    for p in PROCESSES:
        p = dict(p)
        p['species_roles'] = dict(p.get('species_roles', {}))
        if p['requires'] == 'ssu':
            p['requires'] = [{'rule': 'all', 'species': ribo['ssu_proteins']}, {'rule': 'any', 'species': ['R_0069', 'R_0534']}]
            for s in ribo['ssu_proteins']:
                p['species_roles'][s] = '30S ribosomal protein'
            for s in ['R_0069', 'R_0534']:
                p['species_roles'][s] = '16S rRNA (either copy)'
        elif p['requires'] == 'lsu':
            p['requires'] = [{'rule': 'all', 'species': ribo['lsu_proteins']}, {'rule': 'any', 'species': ['R_0068', 'R_0533']},
                             {'rule': 'any', 'species': ['R_0067', 'R_0532']}]
            for s in ribo['lsu_proteins']:
                p['species_roles'][s] = '50S ribosomal protein'
            for s in ['R_0068', 'R_0533']:
                p['species_roles'][s] = '23S rRNA (either copy)'
            for s in ['R_0067', 'R_0532']:
                p['species_roles'][s] = '5S rRNA (either copy)'
        if p['id'] == 'membrane_insertion':
            p['applies_to'] = sorted(g['locus'] for g in genes.values() if g.get('membrane_insertion'))
        out[p['id']] = p
    return out
