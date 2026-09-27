"""Parse input_data/ into the ModelSpec dict.

Every rule here mirrors a line of the simulator, named in the comment beside it; tests/test_consistency.py runs the
simulator's own functions in the container and compares. Values carry a `src` record {file, sheet, row, col} (row = the
spreadsheet row a person would edit, header = row 1) so the browser can say where a number comes from.
"""

import hashlib
import json
import math
import os
import re
import subprocess
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from . import SCHEMA_VERSION
from . import roles as ROLES
from . import names as NAMES

INPUT_FILES = ['initial_concentrations.xlsx', 'kinetic_params.xlsx', 'protein_metabolites.xlsx', 'LargeSubunit.xlsx',
               'oneParamMulder-local_min.json', 'Syn3A_updated.xml', 'syn3A.gb', 'loop_params.txt']

RB_SHEETS = ['Central', 'Nucleotide', 'Lipid', 'Cofactor', 'Transport']     # Rxns_ODE._random_binding_static
H_SPECIES = {'M_h_c', 'M_h_e', 'M_h2o_c', 'M_h2o_e'}                           # dropped by Rxns_ODE.getSpecIDs

# MC_RDME_initialization.initSim: 10 nm lattice, 200 nm radius
LATTICE_SPACING = 10e-9
CYTO_RADIUS_SITES = int(round(2.00e-7 / LATTICE_SPACING))
VOLUME_L = (4 / 3) * np.pi * (CYTO_RADIUS_SITES * LATTICE_SPACING) ** 3 * 1000
NA = 6.022e23
RNAP_SPACING = 400


def mM_to_particles(conc):
    """ImportInitialConditions.mMtoPart"""
    return int(round((conc / 1000) * NA * VOLUME_L))


def _clean(v):
    """numpy/pandas scalar -> JSON value (NaN -> None)."""
    if v is None:
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        return None if math.isnan(v) else float(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    return v


def _numeric(v):
    """(value, stored_as_text): several kinetic_params.xlsx sheets hold numbers as text cells ('804.3384'); the simulator
    passes them through unchanged, the spec keeps the number and flags the cell."""
    if isinstance(v, str):
        try:
            return float(v), True
        except ValueError:
            return v, False
    return v, False


def _src(file, sheet, idx, col):
    return {'file': 'input_data/' + file, 'sheet': sheet, 'row': int(idx) + 2, 'col': col}


def _sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def fingerprint(head):
    """sha256 of every input file; a perturbation records it so a file built against other inputs is flagged."""
    files = {f: _sha(os.path.join(head, 'input_data', f)) for f in INPUT_FILES}
    h = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()[:16]
    return h, files


def _git(head, *args):
    try:
        return subprocess.run(['git', '-C', head] + list(args), capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None


# ------------------------------------------------------------------------------------------------ genome
def load_genome(head):
    """RegionsAndComplexes.mapDNA + createRibosomalProteinMap, without the lattice rotation (kept: original bp coords)."""
    from Bio import SeqIO
    genome = next(SeqIO.parse(os.path.join(head, 'input_data', 'syn3A.gb'), 'gb'))
    genes, ribo_map = {}, {}
    for feat in genome.features:
        q = feat.qualifiers
        if feat.type == 'CDS' and 'protein_id' in q:
            typ = 'protein'
        elif feat.type in ('tRNA', 'rRNA'):
            typ = feat.type
        else:
            continue
        locus = q['locus_tag'][0]
        rna = str(feat.location.extract(genome.seq).transcribe())
        g = {'locus': locus, 'num': locus.split('_')[1], 'type': typ, 'product': str(q['product'][0]),
             'start': int(feat.location.start) + 1, 'end': int(feat.location.end), 'strand': int(feat.strand),
             'rna_len': len(rna), 'src': {'file': 'input_data/syn3A.gb', 'feature': feat.type, 'locus_tag': locus}}
        if 'gene' in q:
            g['symbol'] = q['gene'][0]
        if typ == 'protein':
            g['aa_len'] = len(str(feat.location.extract(genome.seq).transcribe().translate(table=4)).rstrip('*'))
            if 'S ribosomal' in g['product']:
                ribo_map[g['product'].split(' ')[-1]] = locus
        # MC_RDME_initialization.addGeneticInformationReactions: rRNA longer than 2*rnap_spacing -> transcriptionLong
        g['long_transcription'] = typ == 'rRNA' and len(rna) > 2 * RNAP_SPACING
        genes[locus] = g
    return genes, ribo_map, len(genome.seq)


# ------------------------------------------------------------------------------------------------ proteins, RNA
def placed_count(n, localization):
    """ImportInitialConditions.initializeProteins distributes int(n) (membrane) or int(n/2)+int(n/6)+int(n/3)."""
    if localization in ('trans-membrane', 'peripheral membrane'):
        return int(n)
    return int(n / 2) + int(n / 6) + int(n / 3)


def promoter_strength(gene, product, sim_count):
    """ImportInitialConditions.initializePromoterStrengths"""
    if gene['type'] == 'protein':
        if 'S ribosomal' in product:
            return min(765, 500 + sim_count), 'ribosomal: min(765, 500 + count)'
        s = min(765, max(45, sim_count))
        rule = 'min(765, max(45, count))'
        if sim_count < 45:
            rule += ' - floor of 45 applies'
        return s, rule
    if gene['type'] == 'tRNA':
        return 765, 'tRNA: fixed 765'
    return 765 * 3.7, 'rRNA: fixed 765 x 3.7'


def add_expression(head, genes):
    ic = pd.read_excel(os.path.join(head, 'input_data', 'initial_concentrations.xlsx'), sheet_name='Comparative Proteomics')
    for idx, row in ic.iterrows():
        if idx == 0:                                           # row 0 holds citations (initializeProteins skips it)
            continue
        locus = row['Locus Tag']
        g = genes.get(locus)
        if g is None:
            continue
        sim = _clean(row['Sim. Initial Ptn Cnt'])
        g.update({'gene_name': _clean(row['Gene Name']), 'product_proteomics': row['Gene Product'],
                  'essentiality': _clean(row['Essentiality']), 'function': _clean(row['Primary Function']),
                  'localization': _clean(row['Localization']), 'exp_count': _clean(row['Exp. Ptn Cnt']),
                  'sim_count': sim,
                  'src_counts': _src('initial_concentrations.xlsx', 'Comparative Proteomics', idx, 'Sim. Initial Ptn Cnt'),
                  'orthologs': {k: _clean(row[c]) for k, c in [('E. coli', 'Ecoli Ptn Cnt'), ('B. subtilis', 'Bsub Ptn Cnt'),
                                                              ('M. florum', 'M.florum Ptn Cnt')]}})
        n = sim
        if 'S ribosomal' in row['Gene Product']:              # initializeProteins: free ribosomal proteins
            n = max(25, int(sim - 500))
            g['initial_rule'] = 'ribosomal protein: max(25, count - 500) free copies'
        g['initial_protein'] = placed_count(n, g['localization'])
        g['membrane_insertion'] = g['localization'] == 'trans-membrane'   # C_P_ precursor + SecY (initializeProteins)
    for g in genes.values():
        if g['type'] == 'protein' and 'sim_count' not in g:
            g['promoter'], g['promoter_rule'] = None, 'no Comparative Proteomics row: initializePromoterStrengths would fail'
            continue
        strength, rule = promoter_strength(g, g.get('product_proteomics', ''), g.get('sim_count', 0))
        g['promoter'] = strength
        g['promoter_rule'] = rule

    diff = pd.read_excel(os.path.join(head, 'input_data', 'kinetic_params.xlsx'), sheet_name='Diffusion Coefficient')
    for idx, row in diff.iterrows():
        g = genes.get(row['Locus Tag'])
        if g is None:
            continue
        avg = _clean(row['mRNA Avg Count (#)'])
        g['rna_diffusion'] = _clean(row['Diffusion Coeff. (m^2/s)'])
        if g['type'] == 'protein':
            g['mrna_avg'] = avg
            g['mrna_initial_mean'] = 2 * (avg if avg else 0.0001)      # initializeMRNA: Poisson(2 * avg), avg 0 -> 1e-4
            g['src_mrna'] = _src('kinetic_params.xlsx', 'Diffusion Coefficient', idx, 'mRNA Avg Count (#)')
    for g in genes.values():
        if g['type'] == 'tRNA':
            g['initial_rna'] = int(200 / 3) + int(200 / 6) + int(200 / 2)     # initializeTRNA
            g['trna_aa'] = g['product'].split('-')[1].upper()
    return genes


# ------------------------------------------------------------------------------------------------ metabolites
def load_metabolites(head, sbml_names):
    mets = {}
    f = 'initial_concentrations.xlsx'
    ic = pd.read_excel(os.path.join(head, 'input_data', f), sheet_name='Intracellular Metabolites')
    for idx, row in ic.iterrows():
        mid = 'M_' + row['Met ID']
        c = _clean(row['Init Conc (mM)'])
        mets[mid] = {'id': mid, 'kind': 'metabolite', 'name': row['Metabolite name'], 'compartment': 'cytoplasm',
                     'kegg': _clean(row['KEGG ID']), 'initial_mM': c, 'initial_count': mM_to_particles(c),
                     'src': _src(f, 'Intracellular Metabolites', idx, 'Init Conc (mM)')}
    med = pd.read_excel(os.path.join(head, 'input_data', f), sheet_name='Simulation Medium')
    for idx, row in med.iterrows():
        mid = 'M_' + row['Met ID']
        mets[mid] = {'id': mid, 'kind': 'medium', 'name': row['Metabolite name'], 'compartment': 'medium',
                     'kegg': _clean(row['KEGG ID']), 'medium_mM': _clean(row['Conc (mM)']),
                     'note': 'held fixed: the ODE takes medium species as constant parameters (Rxns_ODE)',
                     'src': _src(f, 'Simulation Medium', idx, 'Conc (mM)')}
    for mid, m in mets.items():
        if mid in sbml_names and not m.get('name'):
            m['name'] = sbml_names[mid]
    return mets


# ------------------------------------------------------------------------------------------------ reactions
def rate_law_mm(n_sub, n_prod):
    """Rxns_ODE.Enzymatic, reproduced as a readable string."""
    sn = ' * '.join('(Sub%d/KmSub%d)' % (i, i) for i in range(1, n_sub + 1))
    pn = ' * '.join('(Prod%d/KmProd%d)' % (i, i) for i in range(1, n_prod + 1))
    sd = ' * '.join('(1 + Sub%d/KmSub%d)' % (i, i) for i in range(1, n_sub + 1))
    pd_ = ' * '.join('(1 + Prod%d/KmProd%d)' % (i, i) for i in range(1, n_prod + 1))
    return 'onoff * Enzyme * (kcatF * %s - kcatR * %s) / (%s + %s - 1)' % (sn, pn, sd, pd_)


def _sbml(head):
    import libsbml
    doc = libsbml.readSBMLFromFile(os.path.join(head, 'input_data', 'Syn3A_updated.xml'))
    model = doc.getModel()
    names = [r.name for r in model.getListOfReactions()]
    species_names = {s.getId(): s.getName() for s in model.getListOfSpecies()}
    return model, names, species_names


def _sbml_participants(model, names, rxn_name):
    """Rxns_ODE.getSpecIDs: the SBML reaction whose *name* equals the sheet's reaction name; H+ and water dropped."""
    if rxn_name not in names:
        return None
    r = model.getReaction(names.index(rxn_name))
    subs = [[x.getSpecies(), float(x.getStoichiometry())] for x in r.getListOfReactants() if x.getSpecies() not in H_SPECIES]
    prods = [[x.getSpecies(), float(x.getStoichiometry())] for x in r.getListOfProducts() if x.getSpecies() not in H_SPECIES]
    return {'substrates': subs, 'products': prods, 'sbml_id': r.getId()}


def _enzyme_requires(enz, gpr):
    """Rxns_ODE.getEnzymeConc: one enzyme; 'or' sums the counts; 'and' takes the minimum; 'default' is a fixed 0.001 mM."""
    parts = str(enz).split('-')
    if len(parts) == 1:
        if parts[0] == 'default':
            return [], 'default'
        return [{'rule': 'all', 'species': parts}], 'single'
    return [{'rule': 'all' if gpr == 'and' else 'any', 'species': parts}], gpr


def load_reactions(head, model, sbml_rxn_names):
    rxns = {}
    f = 'kinetic_params.xlsx'
    path = os.path.join(head, 'input_data', f)

    def random_binding(sheet, active):
        df = pd.read_excel(path, sheet_name=sheet)
        for rid, grp in df.groupby('Reaction Name', sort=False):
            r = {'id': rid, 'kind': 'ode_mm', 'layer': 'ODE', 'active': active, 'sheet': sheet,
                 'subsystem': _clean(grp['Subsystem'].iloc[0]) if 'Subsystem' in grp else sheet, 'params': []}
            enz, gpr = None, None
            for idx, row in grp.iterrows():
                pt = row['Parameter Type']
                val, as_text = _numeric(_clean(row['Value']))
                p = {'name': pt, 'value': val, 'unit': _clean(row['Units']), 'src': _src(f, sheet, idx, 'Value')}
                if as_text:
                    p['stored_as_text'] = True
                if _clean(row.get('Related Species')) is not None:
                    p['species'] = row['Related Species']
                if pt == 'Eff Enzyme Count':
                    enz = val
                elif pt == 'GPR rule':
                    gpr = val
                elif pt == 'Substrate Catalytic Rate Constant':
                    p['key'] = 'kcatF'
                elif pt == 'Product Catalytic Rate Constant':
                    p['key'] = 'kcatR'
                elif pt == 'Michaelis Menten Constant':
                    p['key'] = 'Km:' + str(p.get('species'))
                r['params'].append(p)
            r['enzyme_str'] = enz
            r['requires'], r['gpr'] = _enzyme_requires(enz, gpr)
            part = _sbml_participants(model, sbml_rxn_names, rid)
            if part is None:
                r['substrates'], r['products'] = [], []
                r['warning'] = 'no SBML reaction named %r' % rid
            else:
                r['substrates'], r['products'], r['sbml_id'] = part['substrates'], part['products'], part['sbml_id']
                n_sub = int(sum(s for _, s in r['substrates']))
                n_prod = int(sum(s for _, s in r['products']))
                r['rate_law'] = rate_law_mm(n_sub, n_prod)
            kf = next((p['value'] for p in r['params'] if p.get('key') == 'kcatF'), None)
            kr = next((p['value'] for p in r['params'] if p.get('key') == 'kcatR'), None)
            r['reversible'] = bool(kr)
            r['kcatF'], r['kcatR'] = kf, kr
            if not active:
                r['note'] = ('not in the model: defineRxns has the call to defineOtherRandomBindingReactions commented out, '
                             'so this sheet is never read into the ODE')
            rxns[rid] = r

    for sheet in RB_SHEETS:
        random_binding(sheet, True)
    random_binding('Other-Random-Binding', False)

    # Rxns_ODE.defineNonRandomBindingRxns: formula + kinetic law + named parameters; SubN/ProdN rows name species
    sheet = 'Non-Random-Binding Reactions'
    df = pd.read_excel(path, sheet_name=sheet)
    for rid, grp in df.groupby('Reaction Name', sort=False):
        r = {'id': rid, 'kind': 'ode_custom', 'layer': 'ODE', 'active': True, 'sheet': sheet, 'subsystem': 'Non-random binding',
             'params': [], 'substrates': [], 'products': [], 'requires': [], 'gpr': 'none'}
        for idx, row in grp.iterrows():
            pt, val = row['Parameter Type'], _clean(row['Value'])
            if pt == 'Reaction Formula':
                r['formula'] = val
            elif pt == 'Kinetic Law':
                r['rate_law'] = str(val).replace('$', '')
            elif pt.startswith('Sub'):
                r['substrates'].append([val, 1.0])
            elif pt.startswith('Prod'):
                r['products'].append([val, 1.0])
            else:
                num, as_text = _numeric(val)
                r['params'].append({'name': pt, 'key': pt, 'value': num, 'unit': _clean(row['Units']),
                                    'src': _src(f, sheet, idx, 'Value'), **({'stored_as_text': True} if as_text else {})})
        r['reversible'] = bool(next((p['value'] for p in r['params'] if p['name'] == 'kcatR'), 0))
        rxns[rid] = r

    # Rxns_CME.tRNAcharging: synthetase + ATP -> S_atp; + aa -> S_atp_aa; + tRNA -> complex -> charged tRNA + AMP + PPi
    sheet = 'tRNA Charging'
    df = pd.read_excel(path, sheet_name=sheet)
    for rid, grp in df.groupby('Reaction Name', sort=False):
        vals = {row['Parameter Type']: (_clean(row['Value']), idx, _clean(row['Units'])) for idx, row in grp.iterrows()}
        synth, aa = vals['synthetase'][0], vals['amino acid'][0]
        r = {'id': rid, 'kind': 'cme_trna', 'layer': 'CME', 'active': True, 'sheet': sheet, 'subsystem': 'tRNA charging',
             'aa_code': rid[:-3], 'synthetase': synth, 'amino_acid': aa,
             'substrates': [['M_atp_c', 1.0], [aa, 1.0]], 'products': [['M_amp_c', 1.0], ['M_ppi_c', 1.0]],
             'requires': [{'rule': 'all', 'species': [synth]}], 'gpr': 'single',
             'rate_law': 'mass action, 4 steps (k_atp, k_aa, k_tRNA bimolecular; k_cat unimolecular)',
             'params': [{'name': k, 'key': k, 'value': v[0], 'unit': v[2], 'src': _src(f, sheet, v[1], 'Value')}
                        for k, v in vals.items() if k.startswith('k_')]}
        rxns[rid] = r
    return rxns


# ------------------------------------------------------------------------------------------------ protein metabolites
def load_protein_metabolites(head):
    """Rxns_ODE.addProteinMetabolites: the protein's count is split over these ODE species (first = unmodified form)."""
    df = pd.read_excel(os.path.join(head, 'input_data', 'protein_metabolites.xlsx'), sheet_name='protein metabolites')
    out = {}
    for idx, row in df.iterrows():
        out[row['Protein']] = {'forms': row['Metabolite IDs'].split(','),
                               'src': _src('protein_metabolites.xlsx', 'protein metabolites', idx, 'Metabolite IDs')}
    return out


# ------------------------------------------------------------------------------------------------ ribosome assembly
def load_ribosome_assembly(head, ribo_map):
    """Rxns_RDME.addRibosomeBiogenesis: the reduced small-subunit network (<= 19 intermediates) and the large-subunit table."""
    data = json.load(open(os.path.join(head, 'input_data', 'oneParamMulder-local_min.json')))
    max_imts = 19
    sp_names = set(sp['name'] for sp in data['species'] if sp['name'][0] == 'R')
    for err, sps in zip(data['netmin_rmse'], data['netmin_species']):
        if len(sp_names) <= max_imts:
            break
        sp_names.difference_update(set(sps))
    ssu = [r for r in data['reactions'] if r['intermediate'] in sp_names and r['product'] in sp_names]
    rates = {p['rate_id']: p for p in data['parameters'] if p['rate_id'] in set(r['rate_id'] for r in ssu)}
    ssu_steps, ssu_proteins = [], []
    for r in ssu:
        rp = 'S' + r['protein'].split('s')[1]
        locus = ribo_map[rp]
        pid = 'P_' + locus.split('_')[1]
        if pid not in ssu_proteins:
            ssu_proteins.append(pid)
        ssu_steps.append({'protein': pid, 'name': rp, 'intermediate': r['intermediate'], 'product': r['product'],
                          'rate': float(rates[r['rate_id']]['rate']) * 1e6})
    lf = os.path.join(head, 'input_data', 'LargeSubunit.xlsx')
    params = pd.read_excel(lf, sheet_name='parameters')
    lrx = pd.read_excel(lf, sheet_name='reactions')
    lsu_steps, lsu_proteins = [], []
    for idx, row in lrx.iterrows():
        name = row['substrate']
        rate = float(params.loc[params['Protein'] == name]['Rate'].values[0]) * 1e6
        if name == '5S':
            subs = ['R_0067', 'R_0532']
        else:
            locus = ribo_map['L7/L12' if name == 'L7' else name]
            subs = ['P_' + locus.split('_')[1]]
            if subs[0] not in lsu_proteins:
                lsu_proteins.append(subs[0])
        lsu_steps.append({'substrates': subs, 'name': name, 'intermediate': row['intermediate'], 'product': row['product'],
                          'rate': rate, 'src': _src('LargeSubunit.xlsx', 'reactions', idx, 'substrate')})
    return {'ssu_steps': ssu_steps, 'ssu_proteins': ssu_proteins, 'lsu_steps': lsu_steps, 'lsu_proteins': lsu_proteins}


def load_loop_params(head):
    out = {}
    for line in open(os.path.join(head, 'input_data', 'loop_params.txt')):
        line = line.strip()
        if line and not line.startswith('#') and '=' in line:
            k, v = line.split('=', 1)
            out[k.strip()] = v.strip()
    return out


# ------------------------------------------------------------------------------------------------ assembly
def build_spec(head, with_code=True):
    head = os.path.abspath(head)
    genes, ribo_map, genome_len = load_genome(head)
    add_expression(head, genes)
    model, sbml_rxn_names, sbml_species = _sbml(head)
    mets = load_metabolites(head, sbml_species)
    rxns = load_reactions(head, model, sbml_rxn_names)
    ptn_mets = load_protein_metabolites(head)
    ribo = load_ribosome_assembly(head, ribo_map)
    trna_by_aa = {}
    for g in genes.values():
        if g['type'] == 'tRNA':
            trna_by_aa.setdefault(g['trna_aa'], []).append('R_' + g['num'])

    species = {}
    for g in genes.values():
        n = g['num']
        rid = 'R_' + n
        if g['type'] == 'protein':
            pid = 'P_' + n
            species[pid] = {'id': pid, 'kind': 'protein', 'gene': g['locus'], 'name': g.get('product_proteomics') or g['product'],
                            'initial_count': g.get('initial_protein'), 'src': g.get('src_counts'),
                            'compartment': g.get('localization')}
            species[rid] = {'id': rid, 'kind': 'mRNA', 'gene': g['locus'], 'name': 'mRNA of ' + g['locus'],
                            'initial_mean': g.get('mrna_initial_mean'), 'src': g.get('src_mrna')}
        else:
            species[rid] = {'id': rid, 'kind': g['type'], 'gene': g['locus'], 'name': g['product'],
                            'initial_count': g.get('initial_rna', 0)}
    for mid, m in mets.items():
        species[mid] = m
    for pid, pm in ptn_mets.items():
        for i, form in enumerate(pm['forms']):
            sp = species.setdefault(form, {'id': form, 'kind': 'protein_form', 'name': form})
            sp.update({'carrier': pid, 'form_index': i, 'src': pm['src'],
                       'note': ('unmodified form: count = protein %s count minus the other forms' % pid) if i == 0
                       else 'modified form of %s, starts at 0 (initializeProteinMetabolites)' % pid})
            if sp.get('kind') in (None, 'metabolite'):
                sp['kind'] = 'protein_form'
        species['P_' + pid[2:]]['carries'] = pm['forms']
    for cid, info in ROLES.COMPLEXES.items():
        species[cid] = dict(info, id=cid, kind='complex')

    # a reaction that uses a protein's metabolite form cannot run without that protein
    form_owner = {f: pid for pid, pm in ptn_mets.items() for f in pm['forms']}
    for r in rxns.values():
        owners = sorted({form_owner[s] for s, _ in r['substrates'] + r['products'] if s in form_owner})
        for o in owners:
            if not any(o in grp['species'] for grp in r['requires']):
                r['requires'] = r['requires'] + [{'rule': 'all', 'species': [o], 'via': 'protein metabolite'}]
        for s, _ in r['substrates'] + r['products']:
            if s not in species:
                species[s] = {'id': s, 'kind': 'metabolite', 'name': sbml_species.get(s, s),
                              'note': 'appears in a reaction but has no row in Intracellular Metabolites'}
    for r in rxns.values():
        if r['kind'] == 'cme_trna':
            r['trnas'] = trna_by_aa.get(r['aa_code'], [])

    processes = ROLES.build_processes(genes, ribo, trna_by_aa)

    NAMES.check(rxns)
    tc = processes['trna_charging']
    tc['reactions'] = sorted(rid for rid, r in rxns.items() if r['kind'] == 'cme_trna')
    for rid in tc['reactions']:
        r = rxns[rid]
        tc['species_roles'][r['synthetase']] = '%s synthetase (%s)' % (r['aa_code'].title(), rid)
        tc['species_roles'][r['amino_acid']] = 'amino acid charged onto tRNA'
        for t in r.get('trnas', []):
            tc['species_roles'][t] = 'tRNA charged by ' + rid
    tc['species_roles'].update({'M_atp_c': 'ATP consumed per charging', 'M_amp_c': 'AMP released', 'M_ppi_c': 'PPi released'})
    oe = processes['ode_enzymes']
    oe['reactions'] = sorted(rid for rid, r in rxns.items() if r['layer'] == 'ODE' and r['active'])
    for rid in oe['reactions']:
        for grp in rxns[rid]['requires']:
            for sp in grp['species']:
                if sp.startswith('P_'):
                    oe['species_roles'].setdefault(sp, 'enzyme in the ODE')
    for rid, r in rxns.items():
        r['name'] = NAMES.REACTIONS.get(rid, rid)
    used = set()
    for r in rxns.values():
        if r['active']:
            used.update(x for x, _ in r['substrates'] + r['products'])
    for sid, sp in species.items():
        if sp['kind'] == 'medium' and sid not in used:
            sp['unused'] = 'no transport reaction in the model uses it, so changing it does nothing'
        elif sp['kind'] == 'metabolite' and sid not in used and sid not in ROLES.HOOK_METABOLITES:
            sp['unused'] = 'listed in Intracellular Metabolites but no reaction in the model uses it; the count stays fixed'

    # RDME / CME reactions and diffusion, recorded from the simulator's own builders (needs the container's jLM + Bio)
    initial_counts = {sid: sp.get('initial_count', 0) for sid, sp in species.items() if sp.get('kind') == 'metabolite'}
    try:
        import sys
        if head not in sys.path:
            sys.path.insert(0, head)
        from . import rdme as RDME
        rdme = RDME.build(head, genes, ribo_map, VOLUME_L, initial_counts, species)
    except ImportError as e:
        print('WARNING: RDME/CME reactions not recorded (%s); run the build in the 4DWCM container' % e)
        rdme = None
    if rdme:
        for f in list(rdme['families'].values()) + list(rdme['cme_families'].values()):
            if f.get('process') in processes:
                processes[f['process']].setdefault('families', []).append(f['id'])
    for pid in ('trna_charging', 'ode_enzymes'):
        for rid in processes[pid].get('reactions', []):
            rxns[rid]['process'] = pid
    # ribosome assembly intermediates (only named in the RDME) become complex species, so they have pages
    ssu_n, lsu_n = ROLES.SSU, ROLES.LSU
    for f in ((rdme['families'] if rdme else {}) or {}).values():
        if f['family'] not in ('ssu_assembly', 'lsu_assembly', 'ribosome_joining'):
            continue
        for inst in f['instances']:
            for x in inst.get('subs', []) + inst.get('prods', []):
                if x in species or x.startswith(('R_', 'P_')) or x == 'ribosomeP':
                    continue
                small = x.startswith('Rs')
                parts = re.findall(r'5S|[sL]\d+(?:/L12)?', x[1:])
                label = ', '.join(('S' + q[1:]) if q.startswith('s') else q for q in parts)
                species[x] = {'id': x, 'kind': 'complex', 'initial_count': 0, 'assembly': 'SSU' if small else 'LSU',
                              'name': ('30S' if small else '50S') + ' assembly intermediate: ' + ('16S' if small else '23S') + ' rRNA + ' + label,
                              'note': 'intermediate of %s subunit assembly; none at t=0' % ('small' if small else 'large')}
    for cid in (ssu_n, lsu_n):
        if cid in species:
            species[cid]['assembly'] = 'SSU' if cid == ssu_n else 'LSU'

    fp, files = fingerprint(head)
    spec = {'schema': SCHEMA_VERSION,
            'meta': {'head': head, 'built': datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC'),
                     'git_commit': _git(head, 'rev-parse', '--short', 'HEAD'), 'git_branch': _git(head, 'rev-parse', '--abbrev-ref', 'HEAD'),
                     'fingerprint': fp, 'input_sha256': files, 'genome_bp': genome_len,
                     'volume_L': VOLUME_L, 'mM_per_particle': 1000 / (NA * VOLUME_L)},
            'genes': genes, 'species': species, 'reactions': rxns, 'processes': processes,
            'protein_metabolites': ptn_mets, 'ribosome_assembly': ribo, 'loop_params': load_loop_params(head),
            'naming': ROLES.NAMING, 'id_parts': NAMES.ID_PARTS, 'rdme': rdme}
    link_roles(spec)
    if with_code:
        from . import coderefs
        coderefs.attach(spec, head)
    return spec


def link_roles(spec):
    """species -> [{ref, role}] over reactions and processes; role texts say what the species does there."""
    roles = {}
    forms_of = {pid: pm['forms'] for pid, pm in spec.get('protein_metabolites', {}).items()}

    def add(sid, ref, role):
        roles.setdefault(sid, []).append({'ref': ref, 'role': role})

    def arrow(r, forms):
        subs = [x for x, _ in r['substrates'] if x in forms]
        prods = [x for x, _ in r['products'] if x in forms]
        if subs and prods:
            return '%s %s %s' % (' + '.join(subs), '⇌' if r.get('reversible') else '→', ' + '.join(prods))
        return ', '.join(subs + prods)
    for rid, r in spec['reactions'].items():
        for s, st in r['substrates']:
            add(s, 'rxn:' + rid, 'substrate' + (' (x%g)' % st if st != 1 else ''))
        for s, st in r['products']:
            add(s, 'rxn:' + rid, 'product' + (' (x%g)' % st if st != 1 else ''))
        for grp in r['requires']:
            members = grp['species']
            for s in members:
                others = [o for o in members if o != s]
                if grp.get('via'):
                    role = 'as ODE form: ' + arrow(r, forms_of.get(s, []))
                elif r['kind'] == 'cme_trna':
                    role = 'synthetase: charges tRNA-%s' % r['aa_code'].title()
                elif grp['rule'] == 'all' and not others:
                    role = 'sole enzyme'
                elif grp['rule'] == 'all':
                    role = 'subunit, with ' + ', '.join(others) + ' (min count)'
                else:
                    role = 'isozyme, or ' + ', '.join(others) + ' (counts add)'
                add(s, 'rxn:' + rid, role)
        for t in r.get('trnas', []):
            add(t, 'rxn:' + rid, 'tRNA charged with ' + r['aa_code'].title())
    for pid, p in spec['processes'].items():
        for s, role in p.get('species_roles', {}).items():
            if pid == 'ode_enzymes':
                n = sum(1 for r in spec['reactions'].values() if r['active'] and any(s in g['species'] for g in r['requires']))
                role = 'enzyme in %d ODE reaction%s' % (n, '' if n == 1 else 's')
            add(s, 'proc:' + pid, role)
    for fid, f in ((spec.get('rdme') or {}).get('families') or {}).items():
        if f.get('per_gene'):
            continue
        for inst in f['instances']:
            for s in inst.get('subs', []):
                add(s, 'fam:' + fid, 'substrate')
            for s in inst.get('prods', []):
                add(s, 'fam:' + fid, 'product')
    for sid, lst in roles.items():
        if sid in spec['species']:
            seen, uniq = set(), []
            for r in lst:
                k = (r['ref'], r['role'])
                if k not in seen:
                    seen.add(k)
                    uniq.append(r)
            spec['species'][sid]['roles'] = uniq
