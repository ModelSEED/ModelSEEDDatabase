#!/usr/bin/env python
import os, sys
from BiochemPy import Compounds, Reactions,InChIs

with open('../Thermodynamics/Compounds_GroupFormation_eQuilibrator_Comparison.txt') as fh:
    header=1
    Reported_Cpds={'GF':0,'EQ':0,'EQR':0}
    for line in fh.readlines():
        if(header==1):
            header-=1
            continue

        array = line.split('\t')
        if('nan' not in array[1] and '1000000' not in array[1]):
            Reported_Cpds['GF']+=1

        if('nan' not in array[2]):
            Reported_Cpds['EQ']+=1
            (dg,dge)=array[2].split('|')
            if(float(dge)>50):
                Reported_Cpds['EQR']+=1
fh.close()
print("Compounds: ",Reported_Cpds)

with open('../Thermodynamics/Reactions_GroupFormation_eQuilibrator_Comparison.txt') as fh:
    header=1
    Reported_Rxns={'GF':0,'EQ':0,'EQR':0,'GFEQ':0}
    for line in fh.readlines():
        if(header==1):
            header-=1
            continue

        array = line.split('\t')
        if('nan' not in array[1] and '1000000' not in array[1]):
            Reported_Rxns['GF']+=1

        if('nan' not in array[2]):
            Reported_Rxns['EQ']+=1
            if('nan' not in array[1] and '1000000' not in array[1]):
                Reported_Rxns['GFEQ']+=1

            (dg,dge)=array[2].split('|')
            if(float(dge)>100):
                Reported_Rxns['EQR']+=1
fh.close()
print("Reactions: ",Reported_Rxns)
print("================\n")

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

compound_counts=dict()
Cpd_Types = ['All', 'Total (GF)', 'Accepted (GF)', 'Total (EQ)', 'Accepted (EQ)', 'Structured', 'Final']
for type in Cpd_Types:
    compound_counts[type]=0

for cpd in compounds_dict:
    cpd_obj = compounds_dict[cpd]

    compound_counts['All']+=1

    if('GF' in cpd_obj['notes']):
        compound_counts['Total (GF)']+=1

    if('EQ' in cpd_obj['notes']):
        compound_counts['Total (EQ)']+=1

    if('GF' in cpd_obj['notes'] or 'EQ' in cpd_obj['notes']):
        compound_counts['Structured']+=1

    if(cpd_obj['deltag'] == 10000000):
        continue

    if('GF' in cpd_obj['notes'] and 'EQU' not in cpd_obj['notes']):
        compound_counts['Accepted (GF)']+=1
        compound_counts['Final']+=1
    elif('EQU' in cpd_obj['notes']):
        compound_counts['Accepted (EQ)']+=1
        compound_counts['Final']+=1

total_cpds = compound_counts['Structured']
print("Compounds Completeness")
for key in Cpd_Types:
    print(key,compound_counts[key],float(compound_counts[key])/float(total_cpds))
print("================\n")

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

reaction_counts=dict()
Rxn_Types = ['All', 'Total (GF)', 'Complete (GF)', 'Accepted (GF)', 'Total (EQ)', 'Complete (EQ)', 'Accepted (EQ)', 'Structured', 'Final']
for type in Rxn_Types:
    reaction_counts[type]=0

print("Reactions Completeness")
for rxn in reactions_dict:
    rxn_obj = reactions_dict[rxn]

    if(reactions_dict[rxn]['status'] == 'EMPTY'):
        continue

    reaction_counts['All']+=1

    if('GFP' in rxn_obj['notes'] or 'GFC' in rxn_obj['notes']):
        reaction_counts['Total (GF)']+=1

    if('EQP' in rxn_obj['notes'] or 'EQC' in rxn_obj['notes'] or 'Accepted (EQ)' in rxn_obj['notes']):
        reaction_counts['Total (EQ)']+=1

    if('GFP' in rxn_obj['notes'] or 'GFC' in rxn_obj['notes'] or \
           'EQP' in rxn_obj['notes'] or 'EQC' in rxn_obj['notes'] or 'Accepted (EQ)' in rxn_obj['notes']):
        reaction_counts['Structured']+=1

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

total_rxns=reaction_counts['Structured']
for key in Rxn_Types:
    print(key,reaction_counts[key],float(reaction_counts[key])/float(total_rxns))
print("================\n")
