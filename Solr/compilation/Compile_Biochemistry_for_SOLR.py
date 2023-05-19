import os.path
import glob
import json

biochem_root = '../../Biochemistry/'
search_path = os.path.join(biochem_root,"compound_*.json")
global_cpds_list = list()
for compounds_file in sorted(glob.glob(search_path)):
	with open(compounds_file) as json_file_handle:
		cpds_list = json.load(json_file_handle)
		for cpd_obj in cpds_list:
			global_cpds_list.append(cpd_obj)

with open('solr_compounds.json','w') as json_file_handle:
	json_file_handle.write(json.dumps(global_cpds_list))

# "stoichiometry": "-1:cpd00001:0:0:\"H2O\";-1:cpd00012:0:0:\"PPi\";2:cpd00009:0:0:\"Phosphate\";1:cpd00067:0:0:\"H+\""
#                   -1:cpd00001:0:0:\"H2O\";-1:cpd00012:0:0:\"PPi\";2:cpd00009:0:0:\"Phosphate\";1:cpd00067:0:0:\"H+\"
search_path = os.path.join(biochem_root,"reaction_*.json")
global_rxns_list = list()
for reactions_file in sorted(glob.glob(search_path)):
	with open(reactions_file) as json_file_handle:
		rxns_list = json.load(json_file_handle)
		for rxn_obj in rxns_list:

			# compile stoichiometry as string
			stoich_cpd_list = list()
			for cpd_obj in rxn_obj['stoichiometry']:
				cpd_str = ":".join([  str(cpd_obj['coefficient']),
									  cpd_obj['compound'],
									  str(cpd_obj['compartment']),
									  "0",
									  '"'+cpd_obj['name']+'"' ])
				stoich_cpd_list.append(cpd_str)
			stoich_str = ";".join(stoich_cpd_list)

			rxn_obj['stoichiometry'] = stoich_str
			global_rxns_list.append(rxn_obj)

with open('solr_reactions.json','w') as json_file_handle:
	json_file_handle.write(json.dumps(global_rxns_list))


print(len(global_rxns_list))
print(global_rxns_list[0])