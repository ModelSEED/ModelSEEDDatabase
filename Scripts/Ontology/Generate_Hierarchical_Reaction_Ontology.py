#!/usr/bin/env python
from BiochemPy import Reactions
from csv import DictReader
import sys

ONTOLOGY_FILE  = 'Ontology_Donors.tsv'

# The ontology file really needs just two columns with the headers 'from' and 'to'
# Its important to note that the 'from' column must contain the parent (and more generic) metabolite
# and the 'to' column contains the child (and more specific/characterized) metabolite
reader = DictReader(open(ONTOLOGY_FILE), dialect='excel-tab')
compounds_dict = dict()
compounds_ft_dict = dict()
for edge in reader:
    if(edge['to'] not in compounds_dict):
        compounds_dict[edge['to']]=set()
    compounds_dict[edge['to']].add(edge['from'])

# Recursive function for going down through n layers of ontological links
def traverse_ontological_network(node, cpd_dict, neighbor_list):
    for neighbor in cpd_dict[node]:
        if(neighbor not in neighbor_list):
            neighbor_list.add(neighbor)
        if(neighbor in cpd_dict):
            traverse_ontological_network(neighbor, cpd_dict, neighbor_list)

# We iterate through the compounds, and create lists of all their generic parents
# Doesn't have to be their immediate parent.
# We do this as we want to explore all possible links between generic and specific reactions
child_parent_sets = dict()
for child in compounds_dict:
    parents_list=set()
    traverse_ontological_network(child, compounds_dict, parents_list)
    child_parent_sets[child]=parents_list

# load reactions and their codes for matching
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
reactions_codes = reactions_helper.generateCodes(reactions_dict)

(file_stub,suffix) = ONTOLOGY_FILE.rsplit('.', 1)
result_file='.'.join([file_stub,'out'])

outfile = open(result_file,'w')
outfile.write('\t'.join(["Original Reaction","Matched Reaction","Original Description","Matched Description","Swapped Compounds",'\n']))
for rxn in sorted(reactions_dict.keys()):
    if(reactions_dict[rxn]["status"] == "EMPTY"):
        continue

    if(reactions_dict[rxn]["is_obsolete"] == 1):
        continue

    # The use of this parseStoichOnt() function assumes that the biochemistry has been rebuilt and
    # does not include any duplicate reagents or badly formatted floats
    rxn_cpds_dict = reactions_helper.parseStoichOnt(reactions_dict[rxn]["stoichiometry"])

    # We use generateOntologyReactionCodes to generate a list of all possible reactions
    # from a single reaction, that could contain the generic neighbors of the child reagents
    result = reactions_helper.generateOntologyReactionCodes(rxn,rxn_cpds_dict,child_parent_sets)

    # We go through the codes and see what matches
    for new_code in result:
        if(new_code in reactions_codes):
            matched_rxn = sorted(reactions_codes[new_code])[0]
            rxn_desc = reactions_dict[rxn]['definition']
            matched_rxn_desc = reactions_dict[matched_rxn]['definition']

            old_new_array=list()
            for (old,new) in result[new_code]:
                old_new_array.append(old+':'+new)
            old_new = ';'.join(old_new_array)

            outfile.write('\t'.join([rxn,matched_rxn,rxn_desc,matched_rxn_desc,old_new,'\n']))
outfile.close()
