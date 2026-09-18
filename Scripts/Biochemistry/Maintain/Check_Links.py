#!/usr/bin/env python
import os
import sys
sys.path.append('../../../Libs/Python/')
import subprocess
import time
import copy
import re
import json
from collections import OrderedDict
from BiochemPy import Reactions, Compounds

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
cpds_aliases_dict = compounds_helper.loadMSAliases()
cpds_names_dict = compounds_helper.loadNames()
structures_dict = compounds_helper.loadStructures(["InChI","SMILE"],["ModelSEED"])

for cpd in cpds_names_dict:
    if(cpd not in compounds_dict):
        print(cpd+" shouldn't be in names_dict")

for cpd in cpds_aliases_dict:
    if(cpd not in compounds_dict):
        print(cpd+" shouldn't be in aliases_dict")

for cpd in structures_dict:
    if(cpd not in compounds_dict):
        print(cpd+" shouldn't be in structures_dict")

# Load Reactions
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
rxns_aliases_dict = reactions_helper.loadMSAliases()
rxns_names_dict = reactions_helper.loadNames()
rxns_ecs_dict = reactions_helper.loadECs()

for rxn in rxns_names_dict:
    if(rxn not in reactions_dict):
        print(rxn+" shouldn't be in names_dict")

for rxn in rxns_aliases_dict:
    if(rxn not in reactions_dict):
        print(rxn+" shouldn't be in aliases_dict")

for rxn in rxns_ecs_dict:
    if(rxn not in reactions_dict):
        print(rxn+" shouldn't be in ecs_dict")
