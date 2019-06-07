import logging
import cobrakbase
import math
from abc import ABCMeta, ABC, abstractmethod
from equilibrator_api import ComponentContribution, Reaction, Q_
#ccache, Q_, ureg
from equilibrator_cache.exceptions import MissingDissociationConstantsException
from equilibrator_cache.exceptions import ParseException
from sqlalchemy.exc import StatementError
import pandas as pd
from tqdm import tqdm_notebook

logger = logging.getLogger(__name__)

manual_alias = {'cpd11584' : 'mnx:MNXM40496',
                'cpd20863' : 'kegg:C19610',
                'cpd00027' : 'kegg:C00031', 
                'cpd00131' : 'kegg:C00150', 
                'cpd00242' : 'kegg:C00288', 
                'cpd00499' : 'kegg:C01231',
                'cpd00225' : 'kegg:C00261', #C00193 + C00261
                                            #1.0 kegg:C01408 = 2.0 kegg:C00193 -11.7087 kilojoule/mole
                'cpd00077' : 'kegg:C00090',
                'cpd00153' : 'kegg:C00180',
                'cpd00430' : 'kegg:C00548',
                'cpd01629' : 'kegg:C02479',
                'cpd15311' : 'mnx:MNXM31702'
               }

DEFAULT_IONIC_STR = "0.25M"
OVERFLOW = 10000000

def get_kegg_metacyc_energies(kegg_df, meta_df):
    energies = {}
    energies['KEGG'] = {}
    for row_id, d in kegg_df.iterrows():
        cpd_id = d[0]
        #smiles = d[5]
        deltag = d[7]
        deltage = d[8]
        energies['KEGG'][cpd_id] = (deltag, deltage)
        
    energies['MetaCyc'] = {}
    for row_id, d in meta_df.iterrows():
        cpd_id = d[0]
        #smiles = d[5]
        deltag = d[7]
        deltage = d[8]
        energies['MetaCyc'][cpd_id] = (deltag, deltage)
    
    return energies




def get_delta_from_formation(rxn, deltas, deltas_e):
    delta_g_comp = 0
    delta_gerr_comp = 0
    for k in rxn.stoichiometry:
        cpd_deltag = deltas[k]
        if pd.isna(deltas[k]):
            cpd_deltag = OVERFLOW
        if cpd_deltag >= OVERFLOW:
            return OVERFLOW, OVERFLOW
        delta_g_comp += deltas[k] * rxn.stoichiometry[k]
        delta_gerr_comp += math.fabs(deltas_e[k] * rxn.stoichiometry[k])
    return delta_g_comp, delta_gerr_comp / 2

class EquilibratorReactionBuilder:
    
    def __init__(self, modelseed, seed_to_kegg = {}, seed_to_mxn = {}, manual_alias = {}):
        self.seed_to_kegg = seed_to_kegg
        self.seed_to_mxn  = seed_to_mxn
        self.manual_alias = manual_alias
        self.modelseed = modelseed
    
    def get_alias(self, seed_id):
        if seed_id in self.manual_alias:
            return self.manual_alias[seed_id]
        if seed_id in self.seed_to_mxn and len(self.seed_to_mxn[seed_id]) == 1:
            return 'mnx:' + set(self.seed_to_mxn[seed_id]).pop()
        if seed_id in self.seed_to_kegg and len(self.seed_to_kegg[seed_id]) == 1:
            return 'kegg:' + set(self.seed_to_kegg[seed_id]).pop()

        return None
    
    def remap2(self, stoichiometry):
        lhs = {}
        rhs = {}
        for e in stoichiometry:
            alias = self.get_alias(e)
            #print(e, alias)
            if not alias == None:
                if stoichiometry[e] < 0:
                    lhs[alias] = math.fabs(stoichiometry[e])
                elif stoichiometry[e] > 0:
                    rhs[alias] = math.fabs(stoichiometry[e])
                else:
                    print('stoich', e)
                    ok = False
            else:
                #print('!', e, seed_to_mxn[e] if e in seed_to_mxn else 'no alias', 
                #              seed_to_kegg[e] if e in seed_to_kegg else 'no alias')
                1

        return lhs, rhs
    
    def remap3(self, cstoichiometry):
        lhs = {}
        rhs = {}
        for e in cstoichiometry: # e is a tuple (seed_id, cmp)
            seed_id = e[0]
            alias = self.get_alias(seed_id)
            #print(e, alias)
            if not alias == None:
                if cstoichiometry[e] < 0:
                    lhs[alias] = math.fabs(cstoichiometry[e])
                elif cstoichiometry[e] > 0:
                    rhs[alias] = math.fabs(cstoichiometry[e])
                else:
                    print('stoich', e)
                    ok = False
            else:
                #print('!', e, seed_to_mxn[e] if e in seed_to_mxn else 'no alias', 
                #              seed_to_kegg[e] if e in seed_to_kegg else 'no alias')
                1

        return lhs, rhs
    
    def build_equilibrator_equation2(self, cstoichiometry):
        lhs, rhs = self.remap3(cstoichiometry)
        #print(len(stoichiometry), len(lhs), len(rhs))
        if len(cstoichiometry) == len(lhs) + len(rhs):
            eq = to_string(lhs, rhs)
            #print(eq)
            return eq
        else:
            raise Exception("unable to translate compounds")
    
    def build_equilibrator_equation(self, rxn):
        stoichiometry = cobrakbase.modelseed.modelseed.get_stoichiometry(rxn)
        lhs, rhs = self.remap2(stoichiometry)
        #print(len(stoichiometry), len(lhs), len(rhs))
        if len(stoichiometry) == len(lhs) + len(rhs):
            eq = to_string(lhs, rhs)
            #print(eq)
            return eq
        else:
            raise Exception(f'unable to translate compounds')
    
    def build(self, seed_id):
        rxn = self.modelseed.reactions[seed_id]
        eq = self.build_equilibrator_equation(rxn)
        return Reaction.parse_formula(eq)
    

def remap(stoichiometry):
    lhs = {}
    rhs = {}
    for e in stoichiometry:
        if e in seed_to_mxn and len(seed_to_mxn[e]) == 1:
            mnx_id = set(seed_to_mxn[e]).pop()
            if stoichiometry[e] < 0:
                lhs['mnx:' + mnx_id] = math.fabs(stoichiometry[e])
            elif stoichiometry[e] > 0:
                rhs['mnx:' + mnx_id] = math.fabs(stoichiometry[e])
            else:
                print('stoich', e)
                ok = False
        else:
            print('!', e, seed_to_mxn[e] if e in seed_to_mxn else 'no alias')
    return lhs, rhs

def to_string(lhs, rhs):
    return ' + '.join([f'{value} {key}' for key, value in lhs.items()]) + \
           " = " + \
           ' + '.join([f'{value} {key}' for key, value in rhs.items()])

def get_deltag(reaction, cc):
    #rxn_seed_id = rxn['id']
    #eq = build_equilibrator_equation(rxn)
    #print(rxn['id'], eq)
    dG0_prime = None
    dGm_prime = None
    uncertainty = None
    
    dG0_prime, uncertainty = cc.standard_dg_prime(reaction)
    dGm_prime, uncertainty = cc.physiological_dg_prime(reaction)
    return dG0_prime, dGm_prime, uncertainty


def run_for_reactions(rxn_ids, cc, modelseed_dev, units = 'kilocal / mole', abundance = None):
    data_deltas = {}
    seed_to_kegg, seed_to_mxn = get_alias_kegg_mnx(modelseed_dev)
    erxn_builder = EquilibratorReactionBuilder(modelseed_dev, seed_to_kegg, seed_to_mxn, manual_alias)
    for rxn_seed_id in tqdm_notebook(rxn_ids, 'reactions'):
    #for rxn_seed_id in rxn_ids:
        #print(rxn_seed_id)
        rxn = modelseed_dev.get_seed_reaction(rxn_seed_id)
        obsolete = rxn.is_obsolete
        transport = rxn.is_transport
        delta_g = rxn.data['deltag']
        #print(f'{rxn_seed_id}, {obsolete}, {transport}, {delta_g}')
        reaction = None
        balanced = None
        if not obsolete and not transport:
            try:
                reaction = erxn_builder.build(rxn_seed_id)
                balanced = reaction.is_balanced()
                if not abundance == None:
                    for c in reaction.keys():
                        if not c.mnx_id in ['MNXM2', 'MNXM1']:
                            reaction.set_abundance(c, abundance)
                            
                    for c in reaction.keys():
                        print(c.mnx_id, reaction.get_abundance(c))
                dG0_prime, dGm_prime, uncertainty = get_deltag(reaction, cc)
                #eq = build_equilibrator_equation(rxn)
                data_deltas[rxn_seed_id] = (dG0_prime.to(units).magnitude, 
                                            dGm_prime.to(units).magnitude,
                                            uncertainty.to(units).magnitude, 
                                            delta_g, 
                                            reaction.is_balanced(), 
                                            'OK')
                #print(rxn_seed_id, dG0_prime)
            except Exception as e:
                data_deltas[rxn_seed_id] = (None, None, None, None, 
                                            balanced,
                                            e)
        else:
            data_deltas[rxn_seed_id] = (None, None, None, None, 
                                        balanced,
                                        'obsolete/transport')
    return data_deltas
    
    
def run_for_compounds(ids, cc, modelseed_dev, units = 'kilocal / mole'):
    data_deltas = {}
    seed_to_kegg, seed_to_mxn = get_alias_kegg_mnx(modelseed_dev)
    eqb = EquilibratorReactionBuilder(modelseed_dev, seed_to_kegg, seed_to_mxn, manual_alias)
    for cpd_seed_id in tqdm_notebook(ids, 'compounds'):
    #for cpd_seed_id in ids:
        cpd = modelseed_dev.compounds[cpd_seed_id]
        delta_g = cpd['deltag']
        alias = eqb.get_alias(cpd_seed_id)
        equi_delta_g0 = None
        equi_delta_gm = None
        equi_delta_ge = None
        error = 'OK'
        #alias = None
        if not alias == None:
            try:
                reaction = Reaction.parse_formula(' = ' + alias)
                dG0_prime, dGm_prime, uncertainty = get_deltag(reaction, cc)
                #dG0_prime, uncertainty = cc.standard_dg_prime(Reaction.parse_formula(' = ' + alias))
                equi_delta_g0 = dG0_prime.to(units).magnitude
                equi_delta_gm = dGm_prime.to(units).magnitude
                equi_delta_ge = uncertainty.to(units).magnitude
                #print(cpd_seed_id, alias, dG0_prime.to_tuple()[0] / 4.184, delta_g, uncertainty.to_tuple()[0] / 4.184)
            except MissingDissociationConstantsException:
                error = 'MissingDissociationConstantsException'
            except ValueError:
                error = 'ValueError'
            except StatementError:
                error = 'StatementError'
            except ParseException:
                error = 'ParseException'
        else:
            error = 'NoAlias'

        data_deltas[cpd_seed_id] = {
            'seed_delta_g' : delta_g,
            'seed_delta_g_error' : cpd['deltagerr'],
            'alias' : alias,
            'name' : cpd['name'],
            'delta_g' : equi_delta_g0,
            'delta_g_error' : equi_delta_ge,
            'dGm' : equi_delta_gm,
            'error' : error
        }
    return data_deltas

def to_dataframe(data_deltas):
    data = {
        'seed_id' : [],
        'dG0' : [],
        'dGe' : [],
        'dGm' : [],
        'name' : [],
        'alias' : [],
        'error' : []
    }
    for e in data_deltas:
        data['seed_id'].append(e)
        data['dG0'].append(data_deltas[e]['delta_g'])
        data['dGe'].append(data_deltas[e]['delta_g_error'])
        data['dGm'].append(data_deltas[e]['dGm'])
        data['alias'].append(data_deltas[e]['alias'])
        data['name'].append(data_deltas[e]['name'])
        data['error'].append(data_deltas[e]['error'])

    df = pd.DataFrame(data)
    return df

def get_alias_kegg_mnx(ms):
    seed_to_mxn = {}
    seed_to_kegg = {}
    for seed_id in ms.compound_aliases:
        aliases = ms.compound_aliases[seed_id]
        for db in aliases:
            if db == 'metanetx.chemical':
                if not seed_id in seed_to_mxn:
                    seed_to_mxn[seed_id] = set()
                for alias in aliases[db]:
                    seed_to_mxn[seed_id].add(alias)
            elif db == 'KEGG':
                for alias in aliases[db]:
                    if alias.startswith('C') and len(alias) == 6:
                        if not seed_id in seed_to_kegg:
                            seed_to_kegg[seed_id] = set()
                        seed_to_kegg[seed_id].add(alias)
    return seed_to_kegg, seed_to_mxn



def get_basic_rev(dg):
    if dg == OVERFLOW:
        return '?'
    if dg > -2 and dg < 2:
        return '='
    if dg < 0:
        return '>'
    return '<'

def calculate_dg_from_formation(rxn, modelseed_dev):
    compounds = {}
    for cpd_id in rxn.lhs:
        compounds[cpd_id] = modelseed_dev.get_seed_compound(cpd_id)
    for cpd_id in rxn.rhs:
        compounds[cpd_id] = modelseed_dev.get_seed_compound(cpd_id)
    
    dG_comp = 0
    dGe_comp = 0
    for cpd_id in rxn.lhs:
        fdG = compounds[cpd_id].deltag
        fdGe = compounds[cpd_id].data['deltagerr']
        #print(cpd_id, rxn.lhs[cpd_id], fdG, fdGe)
        #print(cpd_id, type(rxn.lhs[cpd_id]), type(fdG), type(fdGe))
        if pd.isna(fdG):
            return OVERFLOW, OVERFLOW
        #print(cpd_id, fdG)
        if fdG == OVERFLOW:
            return OVERFLOW, OVERFLOW
        dG_comp += rxn.lhs[cpd_id] * fdG
        #print(cpd_id, rxn.lhs[cpd_id], fdGe)
        dGe_comp += (rxn.lhs[cpd_id] * fdGe)**2
    for cpd_id in rxn.rhs:
        fdG = compounds[cpd_id].deltag
        fdGe = compounds[cpd_id].data['deltagerr']
        #print(cpd_id, rxn.rhs[cpd_id], fdG, fdGe)
        if pd.isna(fdG):
            return OVERFLOW, OVERFLOW
        #print(cpd_id, fdG)
        if fdG == OVERFLOW:
            return OVERFLOW, OVERFLOW
        dG_comp += rxn.rhs[cpd_id] * fdG
        dGe_comp += (rxn.rhs[cpd_id] * fdGe)**2
    
    dGe_comp = dGe_comp ** 0.5
    
    dG_comp = round(dG_comp, 2)
    dGe_comp = round(dGe_comp, 2)
    
    return dG_comp, dGe_comp

def get_cpd_dg_changes(ms, seed_cpd_dg, databases = ['KEGG', 'MetaCyc']):
    count = 0
    cases = {
        'change' : [],
        'to_null' : [],
        'from_null' : [],
    }
    
    for seed_id in ms.compounds:
        cpd = ms.get_seed_compound(seed_id)
        dG0 = cpd.deltag
        dGe = cpd.data['deltagerr']
        if pd.isna(dG0):
            dG0 = 10000000.0
        match = {}

        for db in databases:
            if cpd.id in ms.compound_aliases and db in ms.compound_aliases[cpd.id]:
                for alias in ms.compound_aliases[cpd.id][db]:
                    if db in seed_cpd_dg and alias in seed_cpd_dg[db]:
                        dG_ = seed_cpd_dg[db][alias]
                        dG_ = (dG_[0], dG_[1])
                        if not dG_ in match:
                            match[dG_] = set()
                        match[dG_].add((alias, db))

        if len(match) > 1:            
            #print(cpd.id, match)
            pass
        elif len(match) == 1:
            p = list(match.keys())[0]
            if not p[0] == dG0 or not p[1] == dGe:
                if dG0 == 10000000:
                    cases['from_null'].append((cpd.id, p))
                elif p[0] == 10000000:
                    cases['to_null'].append((cpd.id, p))
                else:
                    cases['change'].append((cpd.id, p))
                count += 1
                #print(cpd.id, dG0, dGe, match)
        #print(seed_cpd_aliases[cpd.id]['MetaCyc'])
        #break
    for k in cases:
        print(k, len(cases[k]))
        

    print(count)
    return cases

def get_rxn_dg_changes(ms):
    cases = {
        'change' : [],
        'to_null' : [],
        'from_null' : [],
        'reverse' : [],
        'err_change' : [],
    }

    for seed_id in ms.reactions:
        rxn = ms.get_seed_reaction(seed_id)
        if not rxn.is_obsolete and not rxn.is_transport:
            #print(rxn.id)
            dG_comp, dGe_comp = calculate_dg_from_formation(rxn, ms)
            dG_seed = rxn.data['deltag']
            dGe_seed = rxn.data['deltagerr']
            if pd.isna(dG_seed):
                dG_seed = OVERFLOW
            dG_seed = round(dG_seed, 2)
            if not dG_seed == dG_comp or not dGe_seed == dGe_comp:
                dG_change = (seed_id, (dG_seed, dG_comp), (dGe_seed, dGe_comp))
                if dG_seed == OVERFLOW and not dG_comp == OVERFLOW:
                    cases['from_null'].append(dG_change)
                elif not dG_seed == OVERFLOW and dG_comp == OVERFLOW:
                    cases['to_null'].append(dG_change)
                elif dG_seed == dG_comp:
                    cases['err_change'].append(dG_change)
                elif math.fabs(dG_seed) == math.fabs(dG_comp):
                    cases['reverse'].append(dG_change)
                else:
                    cases['change'].append(dG_change)
                #print(seed_id, (dG_seed, dG_comp), (dGe_seed, dGe_comp))

    for k in cases:
        print(k, len(cases[k]))
        
    return cases

def read_deltag_data_from_csv(filename):
    df = pd.read_csv(filename, sep='\t', index_col=0)
    rxn_deltas = {}
    for row_id, d in df.iterrows():
        seed_id = d['seed_id']
        dG0 = OVERFLOW
        dGm = OVERFLOW
        dGe = OVERFLOW
        if 'dG0' in d:
            dG0 = d['dG0']
        if 'dGm' in d:
            dGm = d['dGm']
        if 'dGe' in d:
            dGe = d['dGe']
        if pd.isna(dG0) or dG0 == None:
            dG0 = OVERFLOW
        if pd.isna(dGm) or dGm == None:
            dGm = OVERFLOW
        if pd.isna(dGe) or dGe == None or dGe == '(None,)':
            dGe = OVERFLOW
        dG0 = float(dG0)
        dGm = float(dGm)
        dGe = float(dGe)
        rxn_deltas[seed_id] = (dG0, dGm, dGe)
    return rxn_deltas



def update_compouds_tsv(in_file, out_file, modelseed_dev):
    compouds_tsv = []
    with open(in_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            compouds_tsv.append(l.split('\t'))
        
    for compoud in compouds_tsv:
        seed_id = compoud[0]
        if not seed_id == 'id':
            cpd = modelseed_dev.get_seed_compound(seed_id)
            if not cpd == None:
                if cpd.deltag == None:
                    cpd.deltag = OVERFLOW
                if cpd.data['deltagerr'] == None:
                    cpd.data['deltagerr'] = OVERFLOW
                compoud[12] = str(round(cpd.deltag, 3))
                compoud[13] = str(round(cpd.data['deltagerr'], 3))
            else:
                print('not found', seed_id)
    
    with open(out_file, 'w') as f:
        for compoud in compouds_tsv:
            f.write('\t'.join(compoud))
            
def update_reactions_tsv(in_file, out_file, modelseed_dev):
    reactions_tsv = []
    with open(in_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            reactions_tsv.append(l.split('\t'))
    
    for reaction_data in reactions_tsv:
        seed_id = reaction_data[0]
        if not seed_id == 'id':
            rxn = modelseed_dev.get_seed_reaction(seed_id)
            if not rxn == None:
                if rxn.data['deltag'] == None:
                    rxn.data['deltag']
                if rxn.data['deltagerr'] == None:
                    rxn.data['deltagerr'] = OVERFLOW
                reaction_data[8] = rxn.data['reversibility']
                reaction_data[14] = str(round(rxn.data['deltag'], 2))
                reaction_data[15] = str(round(rxn.data['deltagerr'], 2))
                #compoud[12] = str(round(cpd.deltag, 3))
                #compoud[13] = str(round(cpd.data['deltagerr'], 3))
            else:
                print('not found', seed_id)
                
    with open(out_file, 'w') as f:
        for reaction_data in reactions_tsv:
            f.write('\t'.join(reaction_data))
            
def get_error_data(hist_range, rdata1, rdata2, any_error = True, max_rdata1_error = 20, max_rdata2_error = 20, metabolic_rxns = None):
    data = {
        'seed_id' : [],
        'dG1' : [],
        'dG2' : [],
        'dGe1' : [],
        'dGe2' : [],
        'diff' : [],
    }
    
    error_hist = {}
    rev_hist = {}
    for i in range(hist_range):
        error_hist[i] = set()
        rev_hist[i] = set()
    error_hist['x'] = set()
    rev_hist['x'] = set()
    
    seed_ids = set(rdata1.keys()) & set(rdata2.keys())
    for seed_id in seed_ids:
        dG1  = rdata1[seed_id][0]
        dGe1 = rdata1[seed_id][1]
        dG2  = rdata2[seed_id][0]
        dGe2 = rdata2[seed_id][1]
        
        diff = math.fabs(round(dG1, 2) - round(dG2, 2))
        
        if (metabolic_rxns == None or seed_id in metabolic_rxns) and ((dGe1 <= max_rdata1_error and dGe2 <= max_rdata2_error) or any_error):
            rev1 = rdata1[seed_id][2]
            rev2 = rdata2[seed_id][2]
            data['seed_id'].append(seed_id)
            data['dG1'].append(dG1)
            data['dG2'].append(dG2)
            data['dGe1'].append(dGe1)
            data['dGe2'].append(dGe2)
            data['diff'].append(diff)
            
            added = False
            for i in range(hist_range):
                if diff < i:
                    if diff == 0:
                        i = 0
                    error_hist[i].add(seed_id)
                    rev_hist[i].add((seed_id, rev1, rev2))
                    added = True
                    break
            if not added:
                error_hist['x'].add(seed_id)
                rev_hist['x'].add((seed_id, rev1, rev2))


    df = pd.DataFrame(data)
    df = df.set_index('seed_id')
    
    return error_hist, rev_hist, df

def describe_rev_change(data):
    match = 0
    change_to_irr = 0
    change_to_rev = 0
    change_to_dir = 0
    for p in data:
        rev_old = p[1]
        rev_new = p[2]
        if rev_old == rev_new:
            match += 1
        elif not rev_old == rev_new and rev_old == '=':
            change_to_irr += 1
        elif not rev_old == rev_new and rev_new == '=':
            change_to_rev += 1
        else:
            change_to_dir += 1
    return match, change_to_irr, change_to_rev, change_to_dir

def dg_hist_print(hist_range, error_hist, rev_hist):
    for i in range(hist_range):
        n_m, n_irr, n_rev, n_dir = describe_rev_change(rev_hist[i])
        print(i, len(error_hist[i]), n_m, n_irr, n_rev, n_dir)
    n_m, n_irr, n_rev, n_dir = describe_rev_change(rev_hist['x'])
    print('>%d' % hist_range, len(error_hist['x']), n_m, n_irr, n_rev, n_dir)
    
def get_equilibrator_ok_and_balanced(filename):
    df = pd.read_csv(filename, sep='\t', index_col=0)
    seed_ids = set()
    for row_id, d in df.iterrows():
        seed_id = d['seed_id']
        status = d['status']
        balanced = d['balanced']
        if status == 'OK' and balanced:
            seed_ids.add(seed_id)
    return seed_ids

def fix_null_dg(v):
    if pd.isna(v) or v == None or v =='None':
        return OVERFLOW
    return float(v)

def read_deltag_data_from_csv(filename):
    df = pd.read_csv(filename, sep='\t', index_col=0, low_memory=False)
    rxn_deltas = {}
    for seed_id, d in df.iterrows():
        #print(d)
        dG0 = OVERFLOW
        dGe = OVERFLOW
        reversibility = '?'
        if 'deltag' in d:
            dG0 = d['deltag']
        if 'reversibility' in d:
            reversibility = d['reversibility']
        if 'deltagerr' in d:
            dGe = d['deltagerr']
            
        dG0 = fix_null_dg(dG0)
        dGe = fix_null_dg(dGe)
        
        rxn_deltas[seed_id] = (dG0, dGe, reversibility)
        
    return rxn_deltas

def test_validate_zero(a, b, ids = None):
    both = set(a.keys()) & set(b.keys())
    for rxn_id in both:
        dG1  = a[rxn_id][0]
        dG2  = b[rxn_id][0]
        diff = math.fabs(round(dG1, 2) - round(dG2, 2))
        if diff == 0 and (ids == None or rxn_id in ids):
            if not a[rxn_id][2] == b[rxn_id][2]:
                print(rxn_id, round(dG1, 2), round(dG2, 2), a[rxn_id][2], b[rxn_id][2])
                
                
class ReactionEnergyCalculatorFromCompound:
    
    def __init__(self, compound_energy_calculator):
        self.compound_energy_calculator = compound_energy_calculator
    
    def calculate_reaction_energy(self, stoich):
        dG = 0
        dGe = 0
        missing = set()
        for cpd_id, cmp in stoich:
            
            value = stoich[(cpd_id, cmp)]
            
            energy, error, metadata = self.compound_energy_calculator.calculate_compound_energy(cpd_id)
            
            logger.debug("[%s] %s %s %s", cpd_id, value, energy, error)
            if not energy == None:
                dG += value * energy
                dGe += (value * error)**2
            else:
                missing.add(cpd_id)
            
        dGe = dGe ** 0.5
        if dGe == 0:
            dGe == 2.0
        if len(missing) > 0:
            raise Exception('missing energies for ' + ';'.join(missing))

        return dG, dGe, {}
    
class AbstractEnergyCalculatorFromEquilibrator(metaclass=ABCMeta):
    
    def __init__(self, mapping = {}, calculator = None):
        self.calculator = calculator
        if self.calculator == None:
            self.calculator = ComponentContribution(p_h=Q_(7.0), ionic_strength=Q_("0.25M"), temperature=Q_("298.15K"))
        self.mapping = mapping
        
    def get_deltag(self, equilibrator_reaction, units = 'kilocal / mole'):
        dG0_prime = None
        dGm_prime = None
        uncertainty = None

        dG0_prime, uncertainty = self.calculator.standard_dg_prime(equilibrator_reaction)
        dGm_prime, uncertainty = self.calculator.physiological_dg_prime(equilibrator_reaction)
        dG0_prime = dG0_prime.to(units).magnitude
        dGm_prime = dGm_prime.to(units).magnitude
        uncertainty = uncertainty.to(units).magnitude
        return dG0_prime, dGm_prime, uncertainty
    
class ReactionEnergyCalculatorFromEquilibrator(AbstractEnergyCalculatorFromEquilibrator):
    
    def __init__(self, mapping = {}, calculator = None):
        super().__init__(mapping, calculator)
        self.eq_builder = EquilibratorReactionBuilder(None, {}, {}, mapping)
        pass
    
    def calculate_reaction_energy(self, stoich):
        eq = self.eq_builder.build_equilibrator_equation2(stoich)
        equilibrator_reaction = Reaction.parse_formula(eq)
        logger.debug('Equilibrator reaction: %s', equilibrator_reaction)
        dG0_prime, dGm_prime, uncertainty = self.get_deltag(equilibrator_reaction)
        return dG0_prime, uncertainty, {'dGm' : dGm_prime}

class CompoundEnergyCalculatorFromEquilibrator(AbstractEnergyCalculatorFromEquilibrator):
        
    def calculate_compound_energy(self, seed_id):
        logger.debug(seed_id)
        if seed_id in self.mapping:
            logger.debug('[%s] eQuilibrator alias: %s', seed_id, self.mapping[seed_id])
            equilibrator_reaction = Reaction.parse_formula(' = ' + self.mapping[seed_id])
            dG0_prime, dGm_prime, uncertainty = self.get_deltag(equilibrator_reaction)
            return dG0_prime, uncertainty, {'dGm' : dGm_prime}
        return None, None, {}
    
class CompoundEnergyCalculatorFromFile:
    
    def __init__(self, data = {}):
        self.data = data
        
    def calculate_compound_energy(self, seed_id):
        logger.debug("CompoundEnergyCalculatorFromFile %s", seed_id)
        if seed_id in self.data:
            energy, error = self.data[seed_id]
            return energy, error, {}
        return None, None, {}
    
def get_delta_molanalysis(seed_id, ms, ma, prov, databases = ['KEGG', 'MetaCyc']):
    seed_uinchi = None
    cpd = ms.get_seed_compound(seed_id)
    if seed_id in ms.compound_structures and seed_id in ms.compound_aliases:
        pairs = []
        aliases = ms.compound_aliases[seed_id]
        structures = ms.compound_structures[seed_id]
        if 'InChI' in structures:
            seed_uinchi = structures['InChI']
            for db in databases:
                if db in aliases:
                    for db_id in aliases[db]:
                        if db_id in prov[db].index:
                            l = prov[db].loc[db_id]
                            inchi = l['MASS']
                            if inchi == seed_uinchi:
                                if db in ma and db_id in ma[db]:
                                    dG, dGe = ma[db][db_id]
                                    pairs.append((dG, dGe))
                                
                                #print(db_id, inchi, ma[db][db_id])
                                pass
            #print(modelseed_local.compound_aliases[seed_id])
            
        upairs = set(pairs)
        if len(upairs) == 1:
            gibbs = upairs.pop()
            return gibbs
        elif len(upairs) > 1:
            logger.warning('[%s] Multiple match. Unique Structure: [%s]. Match: %s', seed_id, seed_uinchi, upairs)
            
    return (None, None)