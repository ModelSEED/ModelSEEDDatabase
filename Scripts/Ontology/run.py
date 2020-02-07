#!/usr/bin/python
import pandas as pd
import cobrakbase
from stoich_integration import CStoichiometryHashLibrary
from stoich_integration import multi_hasher
from hierarchical_ontology import HierarchicalOntology

MODELSEED_PATH = '../../'
ONTOLOGY_FILE  = 'ontology_donors.tsv'

modelseed_local = cobrakbase.modelseed.from_local(MODELSEED_PATH)

hlib = CStoichiometryHashLibrary(multi_hasher)
for rxn_id in modelseed_local.reactions:
    rxn = modelseed_local.get_seed_reaction(rxn_id)
    cstoichiometry = rxn.cstoichiometry
    if not rxn.is_obsolete:
        hlib.hash_stoichiometry(rxn_id, cstoichiometry)
        
        
df_ontology = pd.read_csv(ONTOLOGY_FILE, sep='\t')
pairs = set()
for row_id, d in df_ontology.iterrows():
    node_from = d['from']
    node_to = d['to']
    pairs.add((node_from, node_to))
g_cpd = HierarchicalOntology.build_ograph(pairs)
ho = HierarchicalOntology(g_cpd, hlib)

g_rxn = ho.generate_reaction_ontology_from_modelseed(modelseed_local)
print(g_rxn.number_of_nodes(), g_rxn.number_of_edges())

for e in g_rxn.edges:
    rxnA = modelseed_local.get_seed_reaction(e[0])
    rxnB = modelseed_local.get_seed_reaction(e[1])
    print("{}\t{}\t{}\t{}".format(e[0], e[1], rxnA.data['definition'], rxnB.data['definition']))