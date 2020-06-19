#!/usr/bin/env python
import os, sys
from BiochemPy import Compounds, Reactions

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
compound_aliases_dict =  compounds_helper.loadMSAliases()
structures_dict = compounds_helper.loadStructures(["SMILE","InChIKey"],["ModelSEED"])

compound_counts = dict()
compound_structures = dict()
compound_bc_identifiers = dict()
sources = {'BioCyc':{},'Model':{}}

column_keys = ['ModelSEED','KEGG','MetaCyc','KEGG Structures','MetaCyc Structures',
               'BiGG','metanetx.chemical','metanetx.reaction','rhea',
               'KEGG/MetaCyc (Integrated)','KEGG/MetaCyc (Integrated Structure)','KEGG/MetaCyc (Unintegrated)','KEGG/MetaCyc (Unintegrated Structure)',
               'BioCyc (Integrated)','BioCyc (Integrated Identifier)','BioCyc (Unintegrated)',
               'Model (Integrated)', 'Model (Unintegrated)',
               'Orphans']
#NB: Orphans have arisen from several stages of database migration and integration. The most common reason for their occurence is that
#the original parent database has dropped them, and we no longer have the provenance

for key in column_keys:
    compound_counts[key]=0

for cpd in compounds_dict:

    if(compounds_dict[cpd]['is_obsolete'] == 1):
        continue

    compound_counts['ModelSEED']+=1

    if(cpd not in compound_aliases_dict and cpd not in structures_dict):
        #This is possible for a few compounds, but should not be
        continue

    for source in ('KEGG','MetaCyc','BiGG','metanetx.chemical'):
        if(source in compound_aliases_dict[cpd]):
            compound_counts[source]+=1

            if(cpd in structures_dict and (source == 'KEGG' or source == 'MetaCyc')):
                compound_counts[source+' Structures']+=1

    KEGG_MetaCyc=False
    if('KEGG' in compound_aliases_dict[cpd] and 'MetaCyc' in compound_aliases_dict[cpd]):
        KEGG_MetaCyc=True
        compound_counts['KEGG/MetaCyc (Integrated)']+=1
        if(cpd in structures_dict):
            compound_counts['KEGG/MetaCyc (Integrated Structure)']+=1
            compound_structures[cpd]=1
            
    elif('KEGG' in compound_aliases_dict[cpd] or 'MetaCyc' in compound_aliases_dict[cpd]):
        KEGG_MetaCyc=True
        compound_counts['KEGG/MetaCyc (Unintegrated)']+=1

        if(cpd in structures_dict):
            compound_counts['KEGG/MetaCyc (Unintegrated Structure)']+=1
            compound_structures[cpd]=1

    BioCyc=False
    for source in compound_aliases_dict[cpd]:
        if('Cyc' in source and source != "MetaCyc"):
            sources['BioCyc'][source]=1
            BioCyc=True

    if(BioCyc is True):
        BC_Identifier = False
        if('MetaCyc' in compound_aliases_dict[cpd]):
            for mc_alias in compound_aliases_dict[cpd]['MetaCyc']:
                for source in compound_aliases_dict[cpd]:
                    if('Cyc' in source):
                        for bc_alias in compound_aliases_dict[cpd][source]:
                            if(mc_alias == bc_alias):
                                BC_Identifier = True
        
        if(BC_Identifier is True):
            compound_counts['BioCyc (Integrated)']+=1
            compound_counts['BioCyc (Integrated Identifier)']+=1
            compound_bc_identifiers[cpd]=1
        elif(KEGG_MetaCyc is True):
            compound_counts['BioCyc (Integrated)']+=1
        else:
            compound_counts['BioCyc (Unintegrated)']+=1

    Model=False
    for source in compound_aliases_dict[cpd]:
        if('Cyc' not in source and source != 'KEGG' and source != 'metanetx.chemical'):
            Model=True
            sources['Model'][source]=1

    if(Model is True):
        if(KEGG_MetaCyc is True or BioCyc is True):
            compound_counts['Model (Integrated)']+=1
        else:
            compound_counts['Model (Unintegrated)']+=1

    if(KEGG_MetaCyc is False and BioCyc is False and Model is False):
        compound_counts['Orphans']+=1

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
reaction_aliases_dict =  reactions_helper.loadMSAliases()

reaction_counts = dict()
for key in column_keys:
    reaction_counts[key]=0

for rxn in reactions_dict:

    if(reactions_dict[rxn]['status'] == 'EMPTY'):
        continue

    if(reactions_dict[rxn]['is_obsolete'] == 1):
        continue

    reaction_counts['ModelSEED']+=1

    if(rxn not in reaction_aliases_dict):
        #This is possible for a few reactions, but should not be
        continue

    for source in ('KEGG','MetaCyc','BiGG','metanetx.reaction','rhea'):
        if(source in reaction_aliases_dict[rxn]):
            reaction_counts[source]+=1

    KEGG_MetaCyc=False
    if('KEGG' in reaction_aliases_dict[rxn] or 'MetaCyc' in reaction_aliases_dict[rxn]):
        KEGG_MetaCyc=True

    BioCyc=False
    for source in reaction_aliases_dict[rxn]:
        if('Cyc' in source and source != "MetaCyc"):
            BioCyc=True

    Model=False
    for source in reaction_aliases_dict[rxn]:
        if('Cyc' not in source and source != 'KEGG' and 'metanetx' not in source and source != 'rhea'):
            Model=True

    complete_structure=True
    complete_bc_identifier=True
    for entry in reactions_dict[rxn]['compound_ids'].split(';'):
        if(entry not in compound_structures):
            complete_structure=False
        if(entry not in compound_bc_identifiers):
            complete_bc_identifier=False

    #Note: 'if KEGG and MetaCyc' means integrated, 'elif KEGG or MetaCyc' means unintegrated
    if('KEGG' in reaction_aliases_dict[rxn] and 'MetaCyc' in reaction_aliases_dict[rxn]):
        reaction_counts['KEGG/MetaCyc (Integrated)']+=1
        if(complete_structure is True):
            reaction_counts['KEGG/MetaCyc (Integrated Structure)']+=1
    elif('KEGG' in reaction_aliases_dict[rxn] or 'MetaCyc' in reaction_aliases_dict[rxn]):
        reaction_counts['KEGG/MetaCyc (Unintegrated)']+=1
        if(complete_structure is True):
            reaction_counts['KEGG/MetaCyc (Unintegrated Structure)']+=1

    #Note: 'if BioCyc and MetaCyc' means integrated, 'elif BioCyc and not MetaCyc' means unintegrated
    if(BioCyc is True and 'MetaCyc' in reaction_aliases_dict[rxn]):
        reaction_counts['BioCyc (Integrated)']+=1
        if(complete_bc_identifier is True):
            reaction_counts['BioCyc (Integrated Identifier)']+=1
    elif(BioCyc is True and 'MetaCyc' not in reaction_aliases_dict[rxn]):
        reaction_counts['BioCyc (Unintegrated)']+=1

    #Note: 'if Model and BioCyc or KEGG/MetaCyc' means integrated, else unintegrated
    if(Model is True and (BioCyc is True or KEGG_MetaCyc is True)):
        reaction_counts['Model (Integrated)']+=1
    elif(Model is True and BioCyc is not True and KEGG_MetaCyc is not True):
        reaction_counts['Model (Unintegrated)']+=1

    if(KEGG_MetaCyc is False and BioCyc is False and Model is False):
        reaction_counts['Orphans']+=1

print("\n=================")
print("For Table 1 and Table 3\n")

print("->Compound Integration:\n")

print("Structures",len(structures_dict))
for key in column_keys:
    print(key,compound_counts[key])

print("\n->Reaction Integration:\n")
for key in column_keys:
    print(key,reaction_counts[key])
print("\n=================\n")



