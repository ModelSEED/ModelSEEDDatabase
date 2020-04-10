#!/usr/bin/env python
import os,sys
from BiochemPy import Compounds, Reactions

Overridden_Compounds=dict()
header=list()
with open(os.path.dirname(__file__)+'/ACPs_Master_Formula_Charge.txt') as fh:
    for line in fh.readlines():
        line=line.strip()
        array=line.split('\t')

        cpd=array.pop(0)
        Overridden_Compounds[cpd]=1

compounds_helper = Compounds()
structures_dict = compounds_helper.loadStructures(["InChIKey"],["ModelSEED"])
aliases_dict = compounds_helper.loadMSAliases()
compounds_dict = compounds_helper.loadCompounds()

# Load reactions
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

cpds_rxns_dict=dict()
for rxn in reactions_dict:
    for cpd in reactions_dict[rxn]['compound_ids'].split(';'):
        if(cpd not in cpds_rxns_dict):
            cpds_rxns_dict[cpd]=dict()
        cpds_rxns_dict[cpd][rxn]=1

InChIKey_Cpds=dict()
for cpd in sorted (structures_dict):
    if('InChIKey' not in structures_dict[cpd]):
        continue
    for structure in structures_dict[cpd]['InChIKey']:
        if(structure not in InChIKey_Cpds):
            InChIKey_Cpds[structure]=list()
        InChIKey_Cpds[structure].append(cpd)

structure_rxn_count=dict()
for structure in InChIKey_Cpds:
    n_struct_rxns=0
    for cpd in InChIKey_Cpds[structure]:
        n_rxns = 0
        if(cpd in cpds_rxns_dict):
            n_rxns=len(cpds_rxns_dict[cpd])
        n_struct_rxns+=n_rxns
    structure_rxn_count[structure]=n_struct_rxns

with open('Compounds_to_Merge_Spreadsheet.txt','w') as fh:
    for structure_count in sorted(structure_rxn_count.items(), key=lambda item: item[1], reverse=True):
        structure=structure_count[0]
        if(len(InChIKey_Cpds[structure])==1):
            continue

        cells=[str(structure_rxn_count[structure]),structure]
        cells.append(compounds_dict[sorted(InChIKey_Cpds[structure])[0]]['name'])
        for cpd in sorted(InChIKey_Cpds[structure]):

            kegg_string=""
            metacyc_string=""
            if(cpd in aliases_dict and 'KEGG' in aliases_dict[cpd]):
                kegg_string="KEGG:"+"|".join(aliases_dict[cpd]['KEGG'])

            if(cpd in aliases_dict and 'MetaCyc' in aliases_dict[cpd]):
                metacyc_string="MetaCyc:"+"|".join(aliases_dict[cpd]['MetaCyc'])

            cell_string = "; ".join([cpd,kegg_string,metacyc_string])
            cells.append(cell_string)
        fh.write("\t".join(cells)+"\n")

fh.close()
