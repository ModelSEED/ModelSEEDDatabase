#!/usr/bin/env python
import os
import sys
import subprocess
import time
import copy
import re
import json
from collections import OrderedDict
from BiochemPy import Reactions, Compounds

arguments = list(sys.argv)
#Pop filename
arguments = arguments[1:]
if(len(arguments) != 1 or os.path.isfile(arguments[0]) is False):
    print("Error: script must be initiated with the path to the edited json file from Print_Compound_to_Disambiguate.py")
    sys.exit()

File=arguments[0]
Dry_Run=1

##########################################################
#
# Load disambiguation details
#
##########################################################

disambiguated_reactions=list()
with open(File) as jf:
    Disambiguated_Compound = json.load(jf)
disambiguated_reactions = Disambiguated_Compound['reactions']

#Load Reactions
#Using reaction provenance and disambiguation details
#We find the list of reactions that belong to either
#the original compound or the newly disambiguated compound
reactions_helper = Reactions()
reaction_names_dict = reactions_helper.loadNames()
reaction_ecs_dict = reactions_helper.loadECs()

for object in disambiguated_reactions:
    original_names=reaction_names_dict[object['from']['id']]
    disambig_names=reaction_names_dict[object['to']['id']]
    for name in object['names']:
        if(object['names'][name]=="true"):
            if(name not in original_names):
                original_names.append(name)
            if(name in disambig_names):
                disambig_names.remove(name)

        if(object['names'][name]=="false"):
            if(name in original_names):
                original_names.remove(name)
            if(name not in disambig_names):
                disambig_names.append(name)
        
        if(object['names'][name]=="both"):
            if(name not in original_names):
                original_names.append(name)
            if(name not in disambig_names):
                disambig_names.append(name)

    original_ecs=reaction_ecs_dict[object['from']['id']]
    disambig_ecs=reaction_ecs_dict[object['to']['id']]
    for ec in object['ecs']:
        if(object['ecs'][ec]=="true"):
            if(ec not in original_ecs):
                original_ecs.append(ec)
            if(ec in disambig_ecs):
                disambig_ecs.remove(ec)

        if(object['ecs'][ec]=="false"):
            if(ec in original_ecs):
                original_ecs.remove(ec)
            if(ec not in disambig_ecs):
                disambig_ecs.append(ec)
        
        if(object['ecs'][ec]=="both"):
            if(ec not in original_ecs):
                original_ecs.append(ec)
            if(ec not in disambig_ecs):
                disambig_ecs.append(ec)

    reaction_names_dict[object['from']['id']]=original_names
    reaction_names_dict[object['to']['id']]=disambig_names
    reaction_ecs_dict[object['from']['id']]=original_ecs
    reaction_ecs_dict[object['to']['id']]=disambig_ecs

reactions_helper.saveNames(reaction_names_dict)
reactions_helper.saveECs(reaction_ecs_dict)

#Scripts to run afterwards
#./Update_Reaction_Aliases.py
