#!/usr/bin/env python
import os, sys
from BiochemPy import Compounds, Reactions,InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

compound_counts=dict()
for type in ('Total (GF)', 'Accepted (GF)', 'Total (EQ)', 'Accepted (EQ)', 'Total', 'Final'):
    compound_counts[type]=0

for cpd in compounds_dict:
    cpd_obj = compounds_dict[cpd]

    if('GF' in cpd_obj['notes']):
        compound_counts['Total (GF)']+=1

    if('EQ' in cpd_obj['notes']):
        compound_counts['Total (EQ)']+=1

    if('GF' in cpd_obj['notes'] or 'EQ' in cpd_obj['notes']):
        compound_counts['Total']+=1

    if(cpd_obj['deltag'] == 10000000):
        continue

    if('GF' in cpd_obj['notes'] and 'EQU' not in cpd_obj['notes']):
        compound_counts['Accepted (GF)']+=1
        compound_counts['Final']+=1
    elif('EQU' in cpd_obj['notes']):
        compound_counts['Accepted (EQ)']+=1
        compound_counts['Final']+=1

total_cpds = compound_counts['Total']
print("Compounds Completeness")
for key in ('Total', 'Total (GF)', 'Accepted (GF)', 'Total (EQ)', 'Accepted (EQ)', 'Final'):
    print(key,compound_counts[key],float(compound_counts[key])/float(total_cpds))
print("================\n")

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

reaction_counts=dict()
for type in ('Total', 'Total (GF)','Complete (GF)','Accepted (GF)', 'Total (EQ)','Complete (EQ)','Accepted (EQ)','Final'):
    reaction_counts[type]=0

print("Reactions Completeness")
for rxn in reactions_dict:
    rxn_obj = reactions_dict[rxn]

    if(reactions_dict[rxn]['status'] == 'EMPTY'):
        continue

    if('GFP' in rxn_obj['notes'] or 'GFC' in rxn_obj['notes']):
        reaction_counts['Total (GF)']+=1

    if('EQP' in rxn_obj['notes'] or 'EQC' in rxn_obj['notes'] or 'Accepted (EQ)' in rxn_obj['notes']):
        reaction_counts['Total (EQ)']+=1

    if('GFP' in rxn_obj['notes'] or 'GFC' in rxn_obj['notes'] or \
           'EQP' in rxn_obj['notes'] or 'EQC' in rxn_obj['notes'] or 'Accepted (EQ)' in rxn_obj['notes']):
        reaction_counts['Total']+=1

    if('GFC' in rxn_obj['notes']):
        reaction_counts['Complete (GF)']+=1

    if('EQC' in rxn_obj['notes']):
        reaction_counts['Complete (EQ)']+=1

    if(rxn_obj['deltag'] == 10000000):
        continue

    if('EQU' in rxn_obj['notes']):
        reaction_counts['Accepted (EQ)']+=1
        reaction_counts['Final']+=1

    if('GFC' in rxn_obj['notes'] and 'EQU' not in rxn_obj['notes']):
        reaction_counts['Accepted (GF)']+=1
        reaction_counts['Final']+=1

total_rxns=reaction_counts['Total']
for key in ('Total', 'Total (GF)','Complete (GF)','Accepted (GF)', 'Total (EQ)','Complete (EQ)','Accepted (EQ)','Final'):
    print(key,reaction_counts[key],float(reaction_counts[key])/float(total_rxns))
print("================\n")
