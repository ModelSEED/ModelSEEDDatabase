import logging
logger = logging.getLogger(__name__)
import copy
import networkx as nx
from itertools import permutations

def get_permutations(n):
    perm = permutations(range(n), n) 
    return list(perm)

class StoichiometryHashLibrary:
    
    def __init__(self, func_hash):
        self.func_hash = func_hash
        self.all_hashes = {}
        self.rxn_to_hash = {}
        
    def hash_stoichiometry(self, rxn_id, stoich):
        hashes = self.func_hash(stoich)
        self.rxn_to_hash[rxn_id] = hashes
        for h in hashes:
            hash_val = hashes[h]
            if not h in self.all_hashes:
                self.all_hashes[h] = {}
            if not hash_val in self.all_hashes[h]:
                self.all_hashes[h][hash_val] = set()
            self.all_hashes[h][hash_val].add(rxn_id)
            
    def match(self, stoichiometry):
        match = {}
        s_hash = self.func_hash(stoichiometry)
        for hash_type in s_hash:
            hash_val = s_hash[hash_type]
            for hash_type_lib in self.all_hashes:
                if hash_val in self.all_hashes[hash_type_lib]:
                    for rxn_id_match in self.all_hashes[hash_type_lib][hash_val]:
                        if not rxn_id_match in match:
                            match[rxn_id_match] = set()
                        match[rxn_id_match].add((hash_type, hash_type_lib, hash_val))
        return match
    
class CStoichiometryHashLibrary(StoichiometryHashLibrary):
    
    def match(self, cstoichiometry):
        match = {}
        cmps = list(set(map(lambda x : x[1], cstoichiometry)))
        index_tests = get_permutations(len(cmps))
        
        for t in index_tests:
            replace = {}
            for i in range(len(cmps)):
                replace[cmps[i]] = str(t[i])
                
            tstoichiometry = copy.deepcopy(cstoichiometry)
            tstoichiometry = dict(map(lambda x : ((x[0][0], replace[x[0][1]]), x[1]), tstoichiometry.items()))
            
            s_hash = self.func_hash(tstoichiometry)
            
            for hash_type in s_hash:
                hash_val = s_hash[hash_type]
                for hash_type_lib in self.all_hashes:
                    if hash_val in self.all_hashes[hash_type_lib]:
                        for rxn_id_match in self.all_hashes[hash_type_lib][hash_val]:
                            if not rxn_id_match in match:
                                match[rxn_id_match] = set()
                            replace_str = ';'.join(map(lambda x : "{}:{}".format(x[0], x[1]), replace.items()))
                            match[rxn_id_match].add((hash_type, hash_type_lib, hash_val, replace_str))
                
        return match

def multi_hasher(s, ex1 = []):
    s_copy = copy.deepcopy(s)
    #print(s_copy)
    std = hash(frozenset(s_copy.items()))
    for a in s_copy:
        s_copy[a] = -1 * s_copy[a]
    rev = hash(frozenset(s_copy.items()))

    s_copy = copy.deepcopy(s)
    delete = []
    for o in ex1:
        for p in s_copy:
            if o == p[0]:
                delete.append(p)
                
    for p in delete:
        del s_copy[p]
        
    #print(s_copy)
    h_std = hash(frozenset(s_copy.items()))
    for a in s_copy:
        s_copy[a] = -1 * s_copy[a]
    h_rev = hash(frozenset(s_copy.items()))
    
    hashes = {
        'std' : std,
        'rev' : rev,
    }
    
    return hashes    

def single_hasher(s, ex1 = []):
    s_copy = copy.deepcopy(s)
    #print(s_copy)
    std = hash(frozenset(s_copy.items()))
    for a in s_copy:
        s_copy[a] = -1 * s_copy[a]
    rev = hash(frozenset(s_copy.items()))

    s_copy = copy.deepcopy(s)
    delete = []
    for o in ex1:
        for p in s_copy:
            if o == p:
                delete.append(p)
                
    #print(f'delete: {delete}')
    for p in delete:
        del s_copy[p]
        
    #print(s_copy)
    h_std = hash(frozenset(s_copy.items()))
    for a in s_copy:
        s_copy[a] = -1 * s_copy[a]
    h_rev = hash(frozenset(s_copy.items()))
    
    hashes = {
        'std' : std,
        'rev' : rev,
        'no_h_std' : h_std,
        'no_h_rev' : h_rev
    }
    return hashes
        

def filter_hash(ccs):
    clusters = []
    for cc in ccs:
        ids = set()
        for id in cc:
            if not id.endswith('@HASH'):
                ids.add(id)
        if len(ids) > 1:
            clusters.append(ids)
    return clusters

def to_universal(s):
    return replace_keys(s, uid_alias_map)

def replace_keys(s, replace_map):
    s_ = {}
    for i in s:
        if i in replace_map:
            s_[replace_map[i]] = s[i]
        else:
            #print('not found', i)
            s_[i] = s[i]
    return s_

def replace_pkeys(s, replace_map):
    s_ = {}
    for p in s:
        i = p[0]
        if i in replace_map:
            s_[(replace_map[i], p[1])] = s[p]
        else:
            #print('not found', i)
            s_[p] = s[p]
    return s_

def hash_it(rxns, hasher, stoich_f = lambda x : x):
    rxn_to_hash = {}
    all_hashes = {}
    for rxn in rxns:
        rxn_id = '{}@{}'.format(rxn.id, rxn.database)
        hashes = hasher(stoich_f(rxn.cstoichiometry))
        rxn_to_hash[rxn_id] = hashes
        for h in hashes:
            hash_val = hashes[h]
            if not h in all_hashes:
                all_hashes[h] = {}
            if not hash_val in all_hashes[h]:
                all_hashes[h][hash_val] = set()
            all_hashes[h][hash_val].add(rxn_id)
    return rxn_to_hash, all_hashes

def hasher(s, ex1 = []):
    s_copy = copy.deepcopy(s)
    #print(s_copy)
    std = hash(frozenset(s_copy.items()))
    for a in s_copy:
        s_copy[a] = -1 * s_copy[a]
    rev = hash(frozenset(s_copy.items()))

    s_copy = copy.deepcopy(s)
    delete = []
    for o in ex1:
        for p in s_copy:
            #print(f'{o} == {p[0]}')
            if o == p[0]:
                delete.append(p)
                
    #print(f'delete: {delete}')
    for p in delete:
        del s_copy[p]
        
    #print(s_copy)
    h_std = hash(frozenset(s_copy.items()))
    for a in s_copy:
        s_copy[a] = -1 * s_copy[a]
    h_rev = hash(frozenset(s_copy.items()))
    
    hashes = {
        'std' : std,
        'rev' : rev,
        'h_std' : h_std,
        'h_rev' : h_rev
    }
    return hashes

def cluster_reactions(rxns, hasher = hasher, stoich_f = lambda x : x):

    rxn_to_hash, all_hashes = hash_it(rxns, hasher, stoich_f)

    g = nx.Graph()
    for h in all_hashes:
        for h_val in all_hashes[h]:
            #print(h, h_val, all_hashes[h][h_val])
            for i in all_hashes[h][h_val]:
                g.add_edge(str(h_val) + '@HASH', i)
                    #print(prev, i)
            #break

    ccs = nx.algorithms.connected_components(g)
    ccs = [cc for cc in ccs]
    ccs = filter_hash(ccs)
    print('clusters:', len(ccs))
    
    return ccs

class ExternalReference:
    
    def __init__(self, id, namespace):
        self.id = id
        self.namespace = namespace
        
    def __str__(self):
        return "{}@{}".format(self.id, self.namespace)
    
    def __repr__(self):
        return "{}@{}".format(self.id, self.namespace)
    
    def __eq__(self, o):
        return self.id == o.id and self.namespace == o.namespace
    
    def __hash__(self):
        return hash(self.id) * 7 + hash(self.namespace)
    
class MapValidator:
    
    def __init__(self):
        self.mapping = {}
        self.valid_identifiers = {}
    
    def load(self, mapping):
        pass
    
    def add_valid_identifiers(self, id, ns):
        if not ns in self.valid_identifiers:
            self.valid_identifiers[ns] = set()
        self.valid_identifiers[ns].add(id)
    
    def add_annotation(self, e, annotation_id, ns, ns_id):
        if not type(e) == ExternalReference:
            logger.warning('Entity must be type of ExternalReference')
            return None
        if e.namespace in self.valid_identifiers:
            if not e.id in self.valid_identifiers[e.namespace]:
                logger.warning('Invalid identifier for [%s]: %s', e.namespace, e.id)
                return None
        if not e in self.mapping:
            self.mapping[e] = {}
        if not annotation_id in self.mapping[e]:
            self.mapping[e][annotation_id] = {}
        if not ns in self.mapping[e][annotation_id]:
            self.mapping[e][annotation_id][ns] = set()
        self.mapping[e][annotation_id][ns].add(ns_id)
    
    @staticmethod
    def add_map_annotation(a, b):
        #print(a)
        #print(b)
        for ns in b:
            if not ns in a:
                a[ns] = set()
            for o in b[ns]:
                a[ns].add(o)

    def expand(self, annotation_id1, annotation_id2, namespace):
        for e in self.mapping:
            if annotation_id1 in self.mapping[e] and \
               namespace in self.mapping[e][annotation_id1]:
                for i in self.mapping[e][annotation_id1][namespace]:
                    iref = ExternalReference(i, namespace)
                    if iref in self.mapping and \
                       annotation_id2 in self.mapping[iref]:
                        target = self.mapping[e][annotation_id1]
                        transfer = self.mapping[iref][annotation_id2]
                        MapValidator.add_map_annotation(target, transfer)

            #print(e, validator.mapping[e])
            #break
    
    @staticmethod
    def matcher(a, b):
        if a == b:
            return {
                'value' : 'MATCH',
                'data' : set(a)
            }
        return {
            'value' : 'FAIL',
            'data' : (set(a), set(b))
        }
        
    def compare(self, annotation_id1, annotation_id2, namespace, matcher = None):
        if matcher == None:
            matcher = MapValidator.matcher
        result = {}
        for e in self.mapping:
            if annotation_id1 in self.mapping[e] and \
               annotation_id2 in self.mapping[e] and \
               namespace in self.mapping[e][annotation_id1] and \
               namespace in self.mapping[e][annotation_id2]:

                actual = self.mapping[e][annotation_id1][namespace]
                expected = self.mapping[e][annotation_id2][namespace]
                result[str(e)] = matcher(actual, expected)

        return result
    
    
class ModelSEEDMapper:
    
    def __init__(self, ms, s_hash, m_hash, cpd_to_seed):
        self.s_hlib = StoichiometryHashLibrary(s_hash)
        self.m_hlib = CStoichiometryHashLibrary(m_hash)
        self.cpd_to_seed = cpd_to_seed
        for rxn_id in ms.reactions:
            rxn = ms.get_seed_reaction(rxn_id)
            cstoichiometry = rxn.cstoichiometry
            cmps = set(map(lambda x : x[1], cstoichiometry))
            if not rxn.is_obsolete and not rxn.is_transport:
                if len(cmps) == 1:
                    stoichiometry = dict(map(lambda x : (x[0][0], x[1]), cstoichiometry.items()))
                    self.s_hlib.hash_stoichiometry(rxn_id, stoichiometry)
            elif not rxn.is_obsolete:
                if len(cmps) > 1:
                    self.m_hlib.hash_stoichiometry(rxn_id, cstoichiometry)
        
        logger.info('Single Compartment Reactions: %d', len(self.s_hlib.rxn_to_hash))
        logger.info('Multi Compartment Reactions: %d', len(self.m_hlib.rxn_to_hash))

    @staticmethod
    def get_cstoichiometry_model(rxn):
        cstoichiometry = {}
        metabolites = rxn.metabolites
        for m in metabolites:
            cstoichiometry[(m.id, m.compartment)] = metabolites[m]
        return cstoichiometry
    
    def remap_compound(self, id):
        if id in self.cpd_to_seed:
            return self.cpd_to_seed[id]
        return id
    
    def match(self, model_rxn):
        
        cstoichiometry = ModelSEEDMapper.get_cstoichiometry_model(model_rxn)
        cmps = list(set(map(lambda x : x[1], cstoichiometry)))
        match = {}
        
        logger.debug("[%s] %s %s", model_rxn.id, model_rxn, cmps)
        if len(cmps) == 1:
            stoichiometry = dict(map(lambda x : (x[0][0], x[1]), cstoichiometry.items()))
            stoichiometry = dict(map(lambda x : (self.remap_compound(x[0]), x[1]), stoichiometry.items()))
            logger.debug("[%s] %s", model_rxn.id, stoichiometry)
            match = self.s_hlib.match(stoichiometry)
        elif len(cmps) > 1:
            cstoichiometry = dict(map(lambda x : ((self.remap_compound(x[0][0]), x[0][1]), x[1]), cstoichiometry.items()))
            logger.debug("[%s] %s", model_rxn.id, cstoichiometry)
            match = self.m_hlib.match(cstoichiometry)
            
        return match