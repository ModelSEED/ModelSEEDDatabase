#!/usr/bin/env python
import os, sys
from BiochemPy import Compounds, Reactions,InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

print("\n================")
print("For Section: \"Computation of thermodynamic properties of ModelSEED compounds and reaction\"\n")

MS_Complete_Structures=dict()
with open("../../../Biochemistry/Structures/Unique_ModelSEED_Structures.txt") as fh:
    for line in fh.readlines():
        line=line.strip()
        array=line.split('\t')

        if("InChI" in array[5]):
            MS_Complete_Structures[array[5]]=1

MNX_Complete_Structures=dict()
with open("../../../Biochemistry/Structures/MetaNetX/chem_prop.tsv") as fh:
    header=1
    for line in fh.readlines():
        if(line[0] == "#"):
            continue

        line=line.strip()
        array=line.split('\t')

        if('InChI' in array[5]):
            MNX_Complete_Structures[array[5]]=1

Shared_Structures=0
for struct in MS_Complete_Structures:
    if(struct in MNX_Complete_Structures):
        Shared_Structures+=1

print("InChI in ModelSEED: ",len(MS_Complete_Structures.keys()))
print("InChI in MetaNetX: ",len(MNX_Complete_Structures.keys()))
print("Shared InChI: ",Shared_Structures,"\n")

with open('../../Thermodynamics/Compounds_GroupContribution_eQuilibrator_Comparison.txt') as fh:
    header=1
    Reported_Cpds={'GC':0,'EQ':0,'EQR':0,'LoEQE':0,'HiEQE':0}
    for line in fh.readlines():
        if(header==1):
            header-=1
            continue

        array = line.split('\t')

        if(compounds_dict[array[0]]["is_obsolete"] == 1):
            continue

        if('nan' not in array[1] and '1000000' not in array[1]):
            Reported_Cpds['GC']+=1

        if('nan' not in array[2]):
            Reported_Cpds['EQ']+=1
            (dg,dge)=array[2].split('|')
            if(float(dge)>50):
                Reported_Cpds['EQR']+=1
            if(float(dge)<5):
                Reported_Cpds['LoEQE']+=1
            if(float(dge)>100):
                Reported_Cpds['HiEQE']+=1
fh.close()

print("Compounds: ")
print("\tGC: ",Reported_Cpds["GC"])
print("\tEQ: ",Reported_Cpds["EQ"])
print("\tRejected EQ: ",Reported_Cpds["EQR"])

with open('../../Thermodynamics/Reactions_GroupContribution_eQuilibrator_Comparison.txt') as fh:
    header=1
    Reported_Rxns={'GC':0,'EQ':0,'EQR':0,'GCEQ':0,'LoEQE':0,'HiEQE':0}
    for line in fh.readlines():
        if(header==1):
            header-=1
            continue

        array = line.split('\t')

        if(array[0] not in reactions_dict or reactions_dict[array[0]]["is_obsolete"] == 1):
            continue

        if('nan' not in array[1] and '1000000' not in array[1]):
            Reported_Rxns['GC']+=1

        if('nan' not in array[2]):
            Reported_Rxns['EQ']+=1
            if('nan' not in array[1] and '1000000' not in array[1]):
                Reported_Rxns['GCEQ']+=1

            (dg,dge)=array[2].split('|')
            if(float(dge)>100):
                Reported_Rxns['EQR']+=1
            if(float(dge)<5):
                Reported_Rxns['LoEQE']+=1
            if(float(dge)>100):
                Reported_Rxns['HiEQE']+=1
fh.close()
print("Reactions: ")
print("\tGC: ",Reported_Rxns["GC"])
print("\tEQ: ",Reported_Rxns["EQ"])
print("\tShared GC & EQ: ",Reported_Rxns["GCEQ"])
print("\tRejected EQ: ",Reported_Rxns["EQR"],"\n")

print("================")
print("For Section: \"Thermodynamics\"\n")

print("Compounds:")
pct = "{0:.2f}".format(float(Reported_Cpds["LoEQE"])/float(Reported_Cpds["EQ"]))
print("\tLow EQ error: ",Reported_Cpds["LoEQE"],pct)
pct = "{0:.2f}".format(float(Reported_Cpds["HiEQE"])/float(Reported_Cpds["EQ"]))
print("\tHigh EQ error: ",Reported_Cpds["HiEQE"],pct,"\n")
print("Reactions:")
pct = "{0:.2f}".format(float(Reported_Rxns["LoEQE"])/float(Reported_Rxns["EQ"]))
print("\tLow EQ error: ",Reported_Rxns["LoEQE"],pct)
pct = "{0:.2f}".format(float(Reported_Rxns["HiEQE"])/float(Reported_Rxns["EQ"]))
print("\tHigh EQ error: ",Reported_Rxns["HiEQE"],pct,"\n")

print("================")
print("For Table 4\n")

compound_counts=dict()
Cpd_Types = ['All', 'Total (GC)', 'Accepted (GC)', 'Total (EQ)', 'Accepted (EQ)', 'Structured', 'Final']
for type in Cpd_Types:
    compound_counts[type]=0

for cpd in compounds_dict:
    if(compounds_dict[cpd]["is_obsolete"] == 1):
        continue

    cpd_obj = compounds_dict[cpd]

    compound_counts['All']+=1

    if('GC' in cpd_obj['notes']):
        compound_counts['Total (GC)']+=1

    if('EQ' in cpd_obj['notes']):
        compound_counts['Total (EQ)']+=1

    if('GC' in cpd_obj['notes'] or 'EQ' in cpd_obj['notes']):
        compound_counts['Structured']+=1

    if(cpd_obj['deltag'] == 10000000):
        continue

    if('GC' in cpd_obj['notes'] and 'EQU' not in cpd_obj['notes']):
        compound_counts['Accepted (GC)']+=1
        compound_counts['Final']+=1
    elif('EQU' in cpd_obj['notes']):
        compound_counts['Accepted (EQ)']+=1
        compound_counts['Final']+=1

total_cpds = compound_counts['Structured']
print("Compounds Completeness")
for key in Cpd_Types:
    pct = "{0:.2f}".format(float(compound_counts[key])/float(total_cpds))
    print(key,compound_counts[key],pct)

reaction_counts=dict()
Rxn_Types = ['All', 'Total (GC)', 'Complete (GC)', 'Accepted (GC)', 'Total (EQ)', 'Complete (EQ)', 'Accepted (EQ)', 'Structured', 'Final']
for type in Rxn_Types:
    reaction_counts[type]=0

print("\nReactions Completeness")
for rxn in reactions_dict:
    if(reactions_dict[rxn]['status'] == 'EMPTY'):
        continue

    if(reactions_dict[rxn]["is_obsolete"] == 1):
        continue

    rxn_obj = reactions_dict[rxn]

    reaction_counts['All']+=1

    if('GCP' in rxn_obj['notes'] or 'GCC' in rxn_obj['notes']):
        reaction_counts['Total (GC)']+=1

    if('EQP' in rxn_obj['notes'] or 'EQC' in rxn_obj['notes'] or 'Accepted (EQ)' in rxn_obj['notes']):
        reaction_counts['Total (EQ)']+=1

    if('GCP' in rxn_obj['notes'] or 'GCC' in rxn_obj['notes'] or \
           'EQP' in rxn_obj['notes'] or 'EQC' in rxn_obj['notes'] or 'Accepted (EQ)' in rxn_obj['notes']):
        reaction_counts['Structured']+=1

    if('GCC' in rxn_obj['notes']):
        reaction_counts['Complete (GC)']+=1

    if('EQC' in rxn_obj['notes']):
        reaction_counts['Complete (EQ)']+=1

    if(rxn_obj['deltag'] == 10000000):
        continue

    if('EQU' in rxn_obj['notes']):
        reaction_counts['Accepted (EQ)']+=1
        reaction_counts['Final']+=1

    if('GCC' in rxn_obj['notes'] and 'EQU' not in rxn_obj['notes']):
        reaction_counts['Accepted (GC)']+=1
        reaction_counts['Final']+=1

total_rxns=reaction_counts['Structured']
for key in Rxn_Types:
    pct = "{0:.2f}".format(float(reaction_counts[key])/float(total_rxns))
    print(key,reaction_counts[key],pct)
print("\n================\n")
