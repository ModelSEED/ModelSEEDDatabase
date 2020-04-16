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
    original_names=list()
    disambig_names=list()
    for name in object['names']:
        if(object['names'][name]=="true" or object['names'][name]=="both"):
            original_names.append(name)
        if(object['names'][name]=="false" or object['names'][name]=="both"):
            disambig_names.append(name)

    reaction_names_dict[object['from']['id']]=original_names
    reaction_names_dict[object['to']['id']]=disambig_names

reactions_helper.saveNames(reaction_names_dict)
reactions_helper.saveECs(reaction_ecs_dict)

#Scripts to run afterwards
#./Update_Reaction_Aliases.py
