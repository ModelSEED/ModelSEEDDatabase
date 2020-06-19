#!/usr/bin/env python
from BiochemPy import Reactions, Compounds

##########################################
# Fixing compounds
##########################################
compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
for cpd in compounds_dict:
    if(compounds_dict[cpd]['notes']=='null'):
        continue
    new_notes=list()
    for note in compounds_dict[cpd]['notes']:
        if(note == 'GF'):
            new_notes.append('GC')
        else:
            new_notes.append(note)
    compounds_dict[cpd]['notes']=new_notes
compounds_helper.saveCompounds(compounds_dict)

###############################################
# Fixing reactions
###############################################
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
for rxn in reactions_dict:
    if(reactions_dict[rxn]['notes']=='null'):
        continue
    new_notes=list()
    for note in reactions_dict[rxn]['notes']:
        if(note == 'GFP'):
            new_notes.append('GCP')
        elif(note == 'GFC'):
            new_notes.append('GCC')
        else:
            new_notes.append(note)
    reactions_dict[rxn]['notes']=new_notes
reactions_helper.saveReactions(reactions_dict)
