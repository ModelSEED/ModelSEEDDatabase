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

print("Compounds Completeness")
for key in sorted(compound_counts.keys()):
    print(key,compound_counts[key],float(compound_counts[key])/float(len(compounds_dict.keys())))
print("================\n")

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

reaction_counts=dict()
for type in ('GF','EQ','GFC','EQC','EQU','GF & EQ','GFC & EQC','GFC & EQU'):
    reaction_counts[type]=0

print("Reactions Completeness")
for rxn in reactions_dict:
    rxn_obj = reactions_dict[rxn]

    if(rxn_obj['deltag'] == 10000000):
        continue

    if('GFP' in rxn_obj['notes'] or 'GFC' in rxn_obj['notes']):
        reaction_counts['GF']+=1

    if('GFC' in rxn_obj['notes']):
        reaction_counts['GFC']+=1

    if('EQP' in rxn_obj['notes'] or 'EQC' in rxn_obj['notes'] or 'EQU' in rxn_obj['notes']):
        reaction_counts['EQ']+=1

    if('EQC' in rxn_obj['notes'] or 'EQU' in rxn_obj['notes']):
        reaction_counts['EQC']+=1

    if('EQU' in rxn_obj['notes']):
        reaction_counts['EQU']+=1

    if(('GFP' in rxn_obj['notes'] or 'GFC' in rxn_obj['notes']) and \
           ('EQP' in rxn_obj['notes'] or 'EQC' in rxn_obj['notes'] or 'EQU' in rxn_obj['notes'])):
        reaction_counts['GF & EQ']+=1

    if('GFC' in rxn_obj['notes'] and ('EQC' in rxn_obj['notes'] or 'EQU' in rxn_obj['notes'])):
        reaction_counts['GFC & EQC']+=1

    if('GFC' in rxn_obj['notes'] and 'EQU' in rxn_obj['notes']):
        reaction_counts['GFC & EQU']+=1

for key in sorted(reaction_counts.keys()):
    print(key,reaction_counts[key],float(reaction_counts[key])/float(len(reactions_dict.keys())))
print("================\n")
