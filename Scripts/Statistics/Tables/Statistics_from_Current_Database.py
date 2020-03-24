#!/usr/bin/env python
import os, sys
from BiochemPy import Compounds, Reactions,InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

print("Statistics from 2019")
out = open('Growth_Stats.tsv','a')
compound_counts={'Generic':0,'Structure':0}
for cpd in compounds_dict:
    if(compounds_dict[cpd]['smiles'] != ""):
        compound_counts['Structure']+=1
    if('R' in compounds_dict[cpd]['formula']):
        compound_counts['Generic']+=1

print(str(len(compounds_dict.keys()))+" compounds")
for entry in ['Structure','Generic']:
    pct = "{0:.2f}".format(float(compound_counts[entry])/float(len(compounds_dict.keys())))
    print(entry,compound_counts[entry],pct)
out.write('\t'.join(['2019','cpd',str(len(compounds_dict.keys())),str(compound_counts['Structure'])])+'\n')

reaction_counts={'Generic':0,'Complete':0,'Balanced':0}
for rxn in reactions_dict:

    if(reactions_dict[rxn]['status'] == "EMPTY"):
        continue

    complete=True
    generic=False

    for entry in reactions_dict[rxn]['compound_ids'].split(';'):

        if(entry in compounds_dict):

            if(compounds_dict[entry]['smiles'] == ""):
                complete=False

            if('R' in compounds_dict[entry]['formula']):
                generic=True

    if(complete is True):
        reaction_counts['Complete']+=1

    if(generic is True):
        reaction_counts['Generic']+=1

    if('OK' in reactions_dict[rxn]['status']):
        reaction_counts['Balanced']+=1

print(str(len(reactions_dict.keys()))+" reactions")
for entry in ['Complete','Balanced','Generic']:
    pct = "{0:.2f}".format(float(reaction_counts[entry])/float(len(reactions_dict.keys())))
    print(entry,reaction_counts[entry],pct)
out.write('\t'.join(['2019','rxn',str(len(reactions_dict.keys())),str(reaction_counts['Balanced'])])+'\n')
