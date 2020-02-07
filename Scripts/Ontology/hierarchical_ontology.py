import itertools
import copy
import networkx as nx
import pandas as pd

def build_graph_from_dict(dict_data):
    g = nx.DiGraph()
    
    for k in dict_data:
        cpd1 = modelseed_local.get_seed_compound(k)
        for v in dict_data[k]:
            cpd2 = modelseed_local.get_seed_compound(v)
            g.add_edge(k, v)
            #g.add_edge(k + " " + cpd1.name, v + " " + cpd2.name)
    return g

class HierarchicalOntology:
    
    def __init__(self, g, hlib, ontology=None, universal_cpds=None):
        if g==None:
            #old stuff using alias instead of modelseed DO NOT USE THIS always provide graph G
            self.ontology = ontology
            self.g = HierarchicalOntology.build_ograph(ontology)
            self.g_bridge, self.bridges = HierarchicalOntology.generate_bridge_edges(self.g)
            self.t_map = HierarchicalOntology.build_translate_map(self.g_bridge, universal_cpds, 'seed.compound')
        else:
            self.g = g
            self.g_bridge, self.bridges = HierarchicalOntology.generate_bridge_edges(self.g)
            self.t_map = {}
            for e in self.g_bridge.edges:
                if not e[0] in self.t_map:
                    self.t_map[e[0]] = set()
                self.t_map[e[0]].add(e[1])
                
        self.hlib = hlib
        self.child_to_parent = HierarchicalOntology.reverse_child_to_parent(self.t_map)
        self.rxn_g = None
    
    def from_csv(f, sep='\t'):
        df = pd.read_csv(f, sep)
        pairs = set()
        for row_id, d in df.iterrows():
            node_from = d['from']
            node_to = d['to']
            pairs.add((node_from, node_to))
        g = HierarchicalOntology.build_ograph(pairs)
        return HierarchicalOntology(g=g)
        
    def get_compounds():
        return None
    
    @staticmethod
    def build_translate_map(g, universal_cpds, database):
        t_map = {}

        alias_to_id = {}
        for row_id, d in universal_cpds.iterrows():
            alias = d['u_alias']
            cpd_ids = d[database]
            if not pd.isna(cpd_ids):
                alias_to_id[alias] = set(cpd_ids.split('#'))

        for src, dst in g.edges:
            #src = p[0]
            #dst = p[1]
            if src in alias_to_id and dst in alias_to_id:
                for from_id in alias_to_id[src]:
                    for target_id in alias_to_id[dst]:
                        if not from_id in t_map:
                            t_map[from_id] = set()
                        t_map[from_id].add(target_id)

        return t_map
        
    @staticmethod
    def reverse_child_to_parent(t_map):
        result = {}
        for p in t_map:
            for c in t_map[p]:
                if not c in result:
                    result[c] = set()
                result[c].add(p)
        return result

    @staticmethod
    def build_ograph(data):
        g = nx.DiGraph()

        for p in data:
            g.add_edge(p[0], p[1])
        return g
    
    @staticmethod
    def generate_bridge_edges(g):
        bridge = {}
        g_bridge = g.copy()
        for n in g.nodes:
            #print(n)
            t = nx.algorithms.dfs_tree(g, n)
            for k in t.nodes:
                if not k == n and not g_bridge.has_edge(n, k):
                    bridge[(n, k)] = True
                    g_bridge.add_edge(n, k)
        nx.set_edge_attributes(g_bridge, name='bridge', values=bridge)
        return g_bridge, bridge
    
    def get_dot_graph(self):
        dot = HierarchicalOntology.translate_dot(self.g_bridge, self.bridges)
        return dot
    
    def generate_reaction_ontology(self, cstoichiometry, test_single_swaps = True):
            result = {}
            single_rep = []
            replacements = []
            for p in cstoichiometry:
                replace_set = []
                cpd_id = p[0]
                if cpd_id in self.child_to_parent:
                    for other_id in self.child_to_parent[cpd_id]:
                        #print(cpd_id, other_id)
                        replace_set.append((cpd_id, other_id))
                        single_rep.append((cpd_id, other_id))
                if len(replace_set) > 0:
                    replacements.append(replace_set)

            if len(replacements) > 0:
                replacements = list(itertools.product(*replacements))
                if test_single_swaps:
                    for swap in single_rep:
                        replacements.append([swap])

                for swaps in replacements:
                    #print('***', swaps)
                    smap = {x:y for x, y in swaps}
                    #print('smap', smap)
                    stoich_swap = copy.deepcopy(cstoichiometry)
                    stoich_swap = {
                        (x if not x in smap else smap[x], c):y 
                        for (x, c), y in stoich_swap.items()
                    }

                    match = self.hlib.match(stoich_swap)
                    for match_id in match:
                        result[match_id] = smap
            return result
    
    def generate_reaction_ontology_old(self, rxn, hash_f, all_hashes, test_single_swaps = True):
        result = {}
        stoich = rxn.cstoichiometry
        single_rep = []
        replacements = []
        for p in stoich:
            replace_set = []
            cpd_id = p[0]
            if cpd_id in self.child_to_parent:
                for other_id in self.child_to_parent[cpd_id]:
                    #print(cpd_id, other_id)
                    replace_set.append((cpd_id, other_id))
                    single_rep.append((cpd_id, other_id))
            if len(replace_set) > 0:
                replacements.append(replace_set)

        if len(replacements) > 0:
            replacements = list(itertools.product(*replacements))
            if test_single_swaps:
                for swap in single_rep:
                    replacements.append([swap])

            for swaps in replacements:
                #print('***', swaps)
                smap = {x:y for x, y in swaps}
                #print('smap', smap)
                stoich_swap = copy.deepcopy(stoich)
                stoich_swap = {
                    (x if not x in smap else smap[x], c):y 
                    for (x, c), y in stoich_swap.items()
                }
                #hash
                hashes = hash_f(stoich_swap)
                #detect
                for h in hashes:
                    h_val = hashes[h]
                    #print(h, h_val)
                    for h_ in all_hashes:
                        if h_val in all_hashes[h_]:
                            for other_id in all_hashes[h_][h_val]:
                                if not other_id in result:
                                    result[other_id] = smap
        return result
    
    def generate_reaction_ontology_from_modelseed(self, ms):
        matches = {}
        for seed_id in ms.reactions:
            rxn = ms.get_seed_reaction(seed_id)
            match = self.generate_reaction_ontology(rxn.cstoichiometry)
            #match = decorate_match(match, modelseed_local)
            if len(match) > 0:
                #seed_id = decorate_id(seed_id, modelseed_local)
                matches[seed_id] = match
                
        self.rxn_g = nx.DiGraph()
        for rxn_id in matches:
            for parent_id in matches[rxn_id]:
                a = parent_id
                b = rxn_id
                if '@' in a:
                    a = a.split('@')[0]
                if '@' in b:
                    b = b.split('@')[0]
                rxn1 = ms.get_seed_reaction(a)
                rxn2 = ms.get_seed_reaction(b)
                if not rxn1.is_obsolete and not rxn2.is_obsolete:
                    self.rxn_g.add_edge(a, b)
                    
        return self.rxn_g
        
    
    @staticmethod
    def translate_dot(g, bridges = []):
        dot = nx.nx_pydot.to_pydot(g)
        for n in dot.get_nodes():
            #print(n.get_name())
            #n.set_label(f'yay!\n{n}')
            n.set_label(f'{n}')
            n.set_color("blue")
            n.set_style("filled")
            n.set_fillcolor("white")

        for e in dot.get_edges():
            src = e.get_source()
            dst = e.get_destination()
            if src[0] == '"' and src[-1] == '"':
                src = src[1:-1]
            if dst[0] == '"' and dst[-1] == '"':
                dst = dst[1:-1]
            if (src, dst) in bridges:
                #e.set_label('!')
                e.set_style('dashed')
        return dot