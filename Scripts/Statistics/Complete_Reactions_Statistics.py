#!/usr/bin/env python
import os, sys
from BiochemPy import Compounds, Reactions,InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

compound_counts=dict()
for type in ('GF','EQ','GF & EQ','GF & EQU'):
    compound_counts[type]=0

for cpd in compounds_dict:
    cpd_obj = compounds_dict[cpd]

    if(cpd_obj['deltag'] == 10000000):
        continue

    if('GF' in cpd_obj['notes']):
        compound_counts['GF']+=1
    if('EQ' in cpd_obj['notes']):
        compound_counts['EQ']+=1
    if('GF' in cpd_obj['notes'] and 'EQ' in cpd_obj['notes']):
        compound_counts['GF & EQ']+=1
    if('GF' in cpd_obj['notes'] and 'EQU' in cpd_obj['notes']):
        compound_counts['GF & EQU']+=1

for key in sorted(compound_counts.keys()):
    print(key,compound_counts[key])

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()


sys.exit()

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
aliases_dict =  compounds_helper.loadMSAliases()
structures_dict = compounds_helper.loadStructures(["SMILE","InChIKey"],["ModelSEED"])

molanalysis_root=os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/ModelSEED/"

molanalysis_dict=dict()
for db in ['KEGG','MetaCyc']:
    file_name = molanalysis_root+db+'_Charged_MolAnalysis.tbl'
    with open(file_name) as file_handle:
        for line in file_handle.readlines():
            line=line.strip()
            mol = line.split('\t')[0]
            molanalysis_dict[mol]=1

InChI_Cpds=dict()
for cpd in compounds_dict:
    #Assumed to be either KEGG or MetaCyc or Both
    if(cpd in aliases_dict and ('KEGG' in aliases_dict[cpd] or 'MetaCyc' in aliases_dict[cpd])):

        molanalysis_complete=False
        if('KEGG' in aliases_dict[cpd]):
            for kegg in aliases_dict[cpd]['KEGG']:
                if(kegg in molanalysis_dict):
                    molanalysis_complete=True

        if('MetaCyc' in aliases_dict[cpd]):
            for metacyc in aliases_dict[cpd]['MetaCyc']:
                if(metacyc in molanalysis_dict):
                    molanalysis_complete=True

        if(molanalysis_complete==False):
            continue

        if(cpd in structures_dict and "InChIKey" in structures_dict[cpd]):
            InChI_Cpds[cpd]=1

            if(compounds_dict[cpd]['is_obsolete']):
                for link in compounds_dict[cpd]['linked_compound'].split(';'):
                    InChI_Cpds[link]=1

print(len(InChI_Cpds))

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

All_InChI_Rxns=dict()
for rxn in reactions_dict:
    if(reactions_dict[rxn]['is_obsolete']):                
        continue

    if(reactions_dict[rxn]['is_transport']):
        continue

    rxn_cpds_array=reactions_helper.parseStoich(reactions_dict[rxn]["stoichiometry"])

    All_InChI=True
    for rgt in rxn_cpds_array:
        if(rgt['compound'] not in InChI_Cpds):
            All_InChI=False

    if(All_InChI==True):
        All_InChI_Rxns[rxn]=1

print(len(All_InChI_Rxns))
