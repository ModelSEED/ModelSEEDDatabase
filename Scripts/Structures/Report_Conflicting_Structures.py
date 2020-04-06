#!/usr/bin/env python
import os
import sys
import json
from BiochemPy import Compounds, Reactions

# Load Compounds
compounds_helper = Compounds()
aliases_dict = compounds_helper.loadMSAliases()

# Load Reactions
ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

# Load ACPs
Overridden_Fields=dict()
header=list()
with open(os.path.dirname(__file__)+'/ACPs_Master_Formula_Charge.txt') as fh:
    for line in fh.readlines():
        line=line.strip()
        array=line.split('\t')

        cpd=array.pop(0)

        if(len(header)==0):
            header=array
            continue

        if(cpd not in Overridden_Fields):
            Overridden_Fields[cpd]=dict()

        for i in range(len(array)):
            if(array[i] == 'null' or array[i] == 10000000):
                continue
            Overridden_Fields[cpd][header[i]]=array[i]

cpds_rxns_dict=dict()
for rxn in Reactions_Dict:
    for cpd in Reactions_Dict[rxn]['compound_ids'].split(';'):
        if(cpd not in cpds_rxns_dict):
            cpds_rxns_dict[cpd]=dict()
        cpds_rxns_dict[cpd][rxn]=1

# Reactions Compound KEGG MetaCyc Name Formula Charge
cpd_conflicts=dict()
with open('Formula_Conflicts.txt') as file_handle:
    for line in file_handle.readlines():
        line=line.strip()
        array=line.split('\t')

        if(array[0] in Overridden_Fields):
            continue

        cpd_conflicts[array[0]]=len(cpds_rxns_dict[array[0]])
#        print("Formula:",array[0],len(cpds_rxns_dict[array[0]]))

KEGG_URL="https://www.genome.jp/dbget-bin/www_bget?"
MetaCyc_URL="https://biocyc.org/META/NEW-IMAGE?object="
#KEGG MetaCyc Name Formula Charge Assigned to Alias to Use Notes? Completed KEGG link MetaCyc link
with open('Compounds_to_Disambiguate_Spreadsheet.txt','w') as fh:
    for rxn_count in sorted(cpd_conflicts.items(), key=lambda item: item[1], reverse=True):
        cpd = rxn_count[0]
        cells=[str(rxn_count[1]),cpd]

        kegg_list=[]
        if(cpd in aliases_dict and 'KEGG' in aliases_dict[cpd]):
            for kegg in aliases_dict[cpd]['KEGG']:
#                kegg_list.append(KEGG_URL+kegg)
                kegg_list.append(kegg)
        kegg_string=', '.join(kegg_list)
        cells.append(kegg_string)

        metacyc_list=[]
        if(cpd in aliases_dict and 'MetaCyc' in aliases_dict[cpd]):
            for metacyc in aliases_dict[cpd]['MetaCyc']:
#                metacyc_list.append(MetaCyc_URL+metacyc)
                metacyc_list.append(metacyc)
        metacyc_string=', '.join(metacyc_list)
        cells.append(metacyc_string)
        
        cell_string = "; ".join([cpd,kegg_string,metacyc_string])
        fh.write('\t'.join(cells)+"\n")
fh.close()
