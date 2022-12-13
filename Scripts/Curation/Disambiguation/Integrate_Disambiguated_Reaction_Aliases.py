#!/usr/bin/env python
import os
import sys
import subprocess
import time
import copy
import re
import json
from collections import OrderedDict
from BiochemPy import Reactions

arguments = list(sys.argv)
#Pop filename
arguments = arguments[1:]
if(len(arguments) != 1 or os.path.isfile(arguments[0]) is False):
    print("Error: script must be initiated with the path to the edited json file from Print_Reaction_to_Disambiguate.py")
    sys.exit()

File=arguments[0]
Dry_Run=0

##########################################################
#
# Load disambiguation details
#
##########################################################

disambiguated_reaction=None
with open(File) as jf:
    disambiguated_reaction = json.load(jf)
old_rxn_id = disambiguated_reaction['from']['id']

#Requires choice of one reaction, whether new or existing
if(len(disambiguated_reaction['to'])>1):
    print("Error: must pick one exisiting reaction from list:")
    for rxn in disambiguated_reaction['to']:
        print("\t",rxn)
    sys.exit()

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
Aliases_Dict = reactions_helper.loadMSAliases()
Names_Dict = reactions_helper.loadNames()

##########################################################
#
# Create new reaction object
# If the newly disambiguated reaction is new to the biochemistry
# Else we create a whole new reaction
#
##########################################################

new_rxn=None
if(len(disambiguated_reaction['to'])==0 or disambiguated_reaction['to'][0]['id'] is None):

    print("Warning: need to write code to define new reaction!")

    #Create New reaction!
    #Find last identifier and increment
    #last_identifier = list(sorted(reactions_dict))[-1]
    #next_identifier = 'rxn'+str(int(re.sub('^rxn','',last_identifier))+1)

    #new_rxn = OrderedDict({ "id":next_identifier,"name":"null","abbreviation":"null","aliases":"null",
    #                        "formula":"null","mass":"10000000","charge":"0",
    #                        "deltag":"10000000","deltagerr":"10000000","pka":"","pkb":"",
    #                        "inchikey":"","smiles":"",
    #                        "is_cofactor":0,"is_core":0,"is_obsolete":0,
    #                        "abstract_reaction":"null","comprised_of":"null","linked_reaction":"null",
    #                        "source":"","ontology":"class:null|context:null"})

    #new_rxn['formula']=disambiguated_reaction['to'][0]['formula']
    #new_rxn['charge']=disambiguated_reaction['to'][0]['charge']
    #new_rxn['mass']=disambiguated_reaction['to'][0]['mass']

else:

    #Merge with Exisiting reaction!
    new_rxn = reactions_dict[disambiguated_reaction['to'][0]['id']]

reactions_dict[new_rxn['id']]=new_rxn

##########################################################
#
# Handle names and aliases
#
##########################################################

#Names
Old_Names=list()
New_Names=list()
names_dict = disambiguated_reaction['from']['names']
for name in sorted(names_dict):
    if(names_dict[name] is True):
        Old_Names.append(name)
    else:
        if(new_rxn['name'] == "null"):
            new_rxn['name']=name
            new_rxn['abbreviation']=name

        New_Names.append(name)

Names_Dict[old_rxn_id]=Old_Names
if(reactions_dict[old_rxn_id]['name'] not in Old_Names):
    reactions_dict[old_rxn_id]['name'] = sorted(Old_Names)[0]
    reactions_dict[old_rxn_id]['abbreviation'] = sorted(Old_Names)[0]

#Retain names if they exist
if new_rxn['id'] in Names_Dict:
    for name in Names_Dict[new_rxn['id']]:
        #Check if they haven't already been retained
        if(name in Old_Names or name in New_Names):
            continue

        New_Names.append(name)

#Aliases
aliases_dict = disambiguated_reaction['from']['aliases']
Old_Aliases=dict()
New_Aliases=dict()
for alias in aliases_dict:
    for source in aliases_dict[alias]:
        if(aliases_dict[alias][source] is True):
            if(source not in Old_Aliases):
                Old_Aliases[source]=list()
            Old_Aliases[source].append(alias)
        else:
            if(source not in New_Aliases):
                New_Aliases[source]=list()

            if(new_rxn['name']=="null"):
                new_rxn['name']=alias
                new_rxn['abbreviation']=alias

            New_Aliases[source].append(alias)

#Retain aliases if they exist
if new_rxn['id'] in Aliases_Dict:
    for source in Aliases_Dict[new_rxn['id']]:
        for alias in Aliases_Dict[new_rxn['id']][source]:
            #Check if they haven't already been retained
            if(source in Old_Aliases and alias in Old_Aliases[source]):
                continue
            if(source in New_Aliases and alias in New_Aliases[source]):
                continue

            #Warn if source still exists
            if(source in Old_Aliases):
                print("Warning: "+source+" source in new reaction still exists in old reaction ("+"|".join(Old_Aliases[source])+")")

            if(source not in New_Aliases):
                New_Aliases[source]=list()
            New_Aliases[source].append(alias)

Aliases_Dict[old_rxn_id]=Old_Aliases
Aliases_Dict[new_rxn['id']]=New_Aliases

##########################################################
#
# Save names, aliases, and reactions
# These are written to separate files
# If the aliases are updated, then the list of structures
# has to be re-composed
#
##########################################################

if(len(New_Names)>0):
    if(new_rxn['name'] not in New_Names):
        new_rxn['name']=New_Names[0]
        new_rxn['abbreviation']=New_Names[0]
    Names_Dict[new_rxn['id']]=New_Names
    if(Dry_Run==0):
        print("Saving update to names")
        reactions_helper.saveNames(Names_Dict)

if(len(New_Aliases)>0):
    if(Dry_Run==0):
        print("Saving update to aliases")
        reactions_helper.saveAliases(Aliases_Dict)

if(Dry_Run==0):
    print("Saving disambiguation of "+old_rxn_id+" as "+new_rxn['id'])
    reactions_helper.saveReactions(reactions_dict)
