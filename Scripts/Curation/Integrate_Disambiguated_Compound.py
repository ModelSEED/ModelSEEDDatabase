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
Dry_Run=0

##########################################################
#
# Load disambiguation details
#
##########################################################

Disambiguated_Compound=None
with open(File) as jf:
    Disambiguated_Compound = json.load(jf)
Old_Cpd_ID = Disambiguated_Compound['from']['id']

#Requires choice of one compound, whether new or existing
if(len(Disambiguated_Compound['to'])>1):
    print("Error: must pick one exisiting compound from list:")
    for cpd in Disambiguated_Compound['to']:
        print("\t",cpd)
    sys.exit()

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
Aliases_Dict = compounds_helper.loadMSAliases()
Names_Dict = compounds_helper.loadNames()
Structures_Dict = compounds_helper.loadStructures(["InChI","SMILE"],["KEGG","MetaCyc"])

##########################################################
#
# Create new compound object
# If the newly disambiguated compound is new to the biochemistry
# Else we create a whole new compound
#
##########################################################

New_Cpd=None
if(len(Disambiguated_Compound['to'])==0 or Disambiguated_Compound['to'][0]['id'] is None):

    #Create New Compound!
    #Find last identifier and increment
    last_identifier = list(sorted(compounds_dict))[-1]
    next_identifier = 'cpd'+str(int(re.sub('^cpd','',last_identifier))+1)

    New_Cpd = OrderedDict({ "id":next_identifier,"name":"null","abbreviation":"null","aliases":"null",
                            "formula":"null","mass":"10000000","charge":"0",
                            "deltag":"10000000","deltagerr":"10000000","pka":"","pkb":"",
                            "inchikey":"","smiles":"",
                            "is_cofactor":0,"is_core":0,"is_obsolete":0,
                            "abstract_compound":"null","comprised_of":"null","linked_compound":"null",
                            "source":"" })

    New_Cpd['formula']=Disambiguated_Compound['to'][0]['formula']
    New_Cpd['charge']=Disambiguated_Compound['to'][0]['charge']
    New_Cpd['mass']=Disambiguated_Compound['to'][0]['mass']

else:

    #Merge with Exisiting Compound!
    New_Cpd = compounds_dict[Disambiguated_Compound['to'][0]['id']]

compounds_dict[New_Cpd['id']]=New_Cpd

##########################################################
#
# Handle names and aliases
#
##########################################################

#Names
Old_Names=list()
New_Names=list()
names_dict = Disambiguated_Compound['from']['names']
for name in sorted(names_dict):
    if(names_dict[name] is True):
        Old_Names.append(name)
    else:
        if(New_Cpd['name'] == "null"):
            New_Cpd['name']=name
            New_Cpd['abbreviation']=name

        New_Names.append(name)

Names_Dict[Old_Cpd_ID]=Old_Names
if(compounds_dict[Old_Cpd_ID]['name'] not in Old_Names):
    compounds_dict[Old_Cpd_ID]['name'] = sorted(Old_Names)[0]
    compounds_dict[Old_Cpd_ID]['abbreviation'] = sorted(Old_Names)[0]

#Retain names if they exist
if New_Cpd['id'] in Names_Dict:
    for name in Names_Dict[New_Cpd['id']]:
        #Check if they haven't already been retained
        if(name in Old_Names or name in New_Names):
            continue

        New_Names.append(name)

#Aliases
#Structures don't need to changing, as they can be re-listed once aliases are updated
aliases_dict = Disambiguated_Compound['from']['aliases']
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

            if(New_Cpd['name']=="null"):
                New_Cpd['name']=alias
                New_Cpd['abbreviation']=alias

            New_Aliases[source].append(alias)

#Retain aliases if they exist
if New_Cpd['id'] in Aliases_Dict:
    for source in Aliases_Dict[New_Cpd['id']]:
        for alias in Aliases_Dict[New_Cpd['id']][source]:
            #Check if they haven't already been retained
            if(source in Old_Aliases and alias in Old_Aliases[source]):
                continue
            if(source in New_Aliases and alias in New_Aliases[source]):
                continue

            #Warn if source still exists
            if(source in Old_Aliases):
                print("Warning: Alias "+alias+" for "+source+" is unresolved")

            if(source not in New_Aliases):
                New_Aliases[source]=list()
            New_Aliases[source].append(alias)

Aliases_Dict[Old_Cpd_ID]=Old_Aliases
Aliases_Dict[New_Cpd['id']]=New_Aliases

##########################################################
#
# Save names, aliases, and compounds
# These are written to separate files
# If the aliases are updated, then the list of structures
# has to be re-composed
#
##########################################################

if(len(New_Names)>0):
    if(New_Cpd['name'] not in New_Names):
        New_Cpd['name']=New_Names[0]
        New_Cpd['abbreviation']=New_Names[0]
    Names_Dict[New_Cpd['id']]=New_Names
    if(Dry_Run==0):
        print("Saving update to names")
        compounds_helper.saveNames(Names_Dict)

if(len(New_Aliases)>0):
    if(Dry_Run==0):
        print("Saving update to aliases")
        compounds_helper.saveAliases(Aliases_Dict)

if(Dry_Run==0):
    print("Saving disambiguation of "+Old_Cpd_ID+" as "+New_Cpd['id'])
    compounds_helper.saveCompounds(compounds_dict)

##########################################################
#
# Find reactions to disambiguate
#
##########################################################

#Reactions
#Read all Provenance
Provenance_Root=os.path.dirname(__file__)+"/../../Biochemistry/Aliases/Provenance/"
Prov_Rxns=dict()
with open(Provenance_Root+"Original_Reaction_Compound_Links.tsv") as pfh:
    for line in pfh.readlines():
        line=line.strip()
        (model,rxn,cpds)=line.split("\t")

        Found_Cpd=""
        for cpd in cpds.split("|"):
            if(cpd in aliases_dict and model in aliases_dict[cpd]):
                Found_Cpd=cpd

        if(Found_Cpd == ""):
            continue

        Is_Old=False
        if(Found_Cpd != "" and model in Old_Aliases and Found_Cpd in Old_Aliases[model]):
            Is_Old=True

        if(model not in Prov_Rxns):
            Prov_Rxns[model]=dict()
        Prov_Rxns[model][rxn]=Is_Old
        
#Load Reactions
#Using reaction provenance and disambiguation details
#We find the list of reactions that belong to either
#the original compound or the newly disambiguated compound
reactions_helper = Reactions()
reaction_aliases_dict = reactions_helper.loadMSAliases()
reaction_names_dict = reactions_helper.loadNames()
reaction_ecs_dict = reactions_helper.loadECs()
reactions_dict = reactions_helper.loadReactions()

#Three types of reactions:
# 1) Retain original/old reaction
keep_reactions=dict()
# 2) Retain but change compound in stoichiometry
update_reactions=dict()
# 3) Split into new reaction with updated stoichiometry
disambiguate_reactions=dict()
for rxn in reaction_aliases_dict:
    Is_Old_Dict=dict()
    Source_Alias=dict()

    for source in reaction_aliases_dict[rxn]:
        if(source not in Prov_Rxns):
            continue

        for alias in reaction_aliases_dict[rxn][source]:
            if(alias not in Prov_Rxns[source]):
                continue

            #This is key, and it happens because we're using external aliases
            #It's entirely possible that:
            #    (1) the reaction alias is associated with multiple distinct ModelSEED reactions
            #    (2) the compound alias, for a given source, was not associated with the compound in question
            if(Old_Cpd_ID not in reactions_dict[rxn]['stoichiometry']):
                continue

            Is_Old_Dict[Prov_Rxns[source][alias]]=1

    if(len(Is_Old_Dict)==1):
        if(list(Is_Old_Dict)[0] is True):
            #This means the reaction is not changed
            keep_reactions[rxn]=1
        else:
            #This means that the reaction's stoichiometry must be updated (in place)
            update_reactions[rxn]=1
    elif(len(Is_Old_Dict)==2):
        #This means that the reaction has previously incorrectly merged reactions
        #And the reaction needs to be split into two new reactions
        disambiguate_reactions[rxn]=dict()

if(len(update_reactions)==0 and len(disambiguate_reactions)==0):
    print("No reactions need updating")
    sys.exit()

##########################################################
#
# Handle reactions to update
#
##########################################################

print("Updating "+str(len(update_reactions))+" reactions")
for rxn in sorted(update_reactions):

    #Parse old stoichiometry into array
    old_stoichiometry=reactions_dict[rxn]["stoichiometry"]
    rxn_cpds_array=reactions_helper.parseStoich(old_stoichiometry)

    #Adjust for new compound
    reactions_helper.replaceCompound(rxn_cpds_array,Old_Cpd_ID,New_Cpd['id'])

    #Rebuild with old compound
    new_stoichiometry = reactions_helper.buildStoich(rxn_cpds_array)
    reactions_helper.rebuildReaction(reactions_dict[rxn],new_stoichiometry)

##########################################################
#
# Handle reactions to disambiguate
# Names and ECs must be disambiguated further in:
# 
##########################################################

#Find last identifier
last_identifier = list(sorted(reactions_dict))[-1]
identifier_count = int(re.sub('^rxn','',last_identifier))

#Generate codes for matching
reactions_codes = reactions_helper.generateCodes(reactions_dict)

print("Disambiguating "+str(len(disambiguate_reactions))+" reactions")
new_reaction_count=dict()
metadata=Disambiguated_Compound['metadata']
disambiguation_object = {'metadata':metadata,
                         'reactions':[],
                         'method':"disambiguation"}

for original_rxn in sorted(disambiguate_reactions):

    #Parse old stoichiometry into array
    old_stoichiometry=reactions_dict[original_rxn]["stoichiometry"]
    rxn_cpds_array=reactions_helper.parseStoich(old_stoichiometry)
    old_rxn_code = reactions_helper.generateCode(rxn_cpds_array)

    #Adjust for new compound and retrieve code for checking against database
    reactions_helper.replaceCompound(rxn_cpds_array,Old_Cpd_ID,New_Cpd['id'])
    rxn_code = reactions_helper.generateCode(rxn_cpds_array)

    disambig_rxn = None
    already_in_database = False
    if(rxn_code in reactions_codes):
        disambig_rxn = sorted(list(reactions_codes[rxn_code]))[0]
        already_in_database = True

    #If there is a match, we assume nothing needs changing
    #Other than aliases. Move aliases from 'new' to 'matched'
    if(disambig_rxn is None):

        #Increment identifier
        identifier_count+=1
        next_identifier = 'rxn'+str(identifier_count)

        #Copy reaction
        new_rxn = copy.deepcopy(reactions_dict[original_rxn])
        new_rxn['id']=next_identifier
        disambig_rxn = new_rxn['id']

        reactions_dict[new_rxn['id']]=new_rxn
        new_reaction_count[new_rxn['id']]=1

        #Rebuild with old compound
        new_stoichiometry = reactions_helper.buildStoich(rxn_cpds_array)
        reactions_helper.rebuildReaction(new_rxn,new_stoichiometry)

        #Finally, because several new reactions may share equations
        if(rxn_code not in reactions_codes):
            reactions_codes[rxn_code]=dict()
        reactions_codes[rxn_code][new_rxn['id']]=1

    #Handle Aliases
    keep_reaction_aliases = dict()
    new_reaction_aliases = dict()
    for source in reaction_aliases_dict[original_rxn]:
        for alias in reaction_aliases_dict[original_rxn][source]:
            if(source in Prov_Rxns and alias in Prov_Rxns[source]):
                #Need to find and handle aliases that are in loaded provenance
                if(Prov_Rxns[source][alias] is True):
                    if(source not in keep_reaction_aliases):
                        keep_reaction_aliases[source]=list()
                    keep_reaction_aliases[source].append(alias)
                else:
                    if(source not in new_reaction_aliases):
                        new_reaction_aliases[source]=list()
                    new_reaction_aliases[source].append(alias)

    #Need to find and retain aliases that are not in loaded provenance
    #This occurs for one or two reasons:
    #    (i) original provenance is missing
    #    (ii) source is an abstraction (i.e. BiGG)
    #Here, we check to see if the alias itself has been moved, regardless of source
    for original_source in reaction_aliases_dict[original_rxn]:
        for original_alias in reaction_aliases_dict[original_rxn][original_source]:
            if(original_source not in Prov_Rxns or original_alias not in Prov_Rxns[original_source]):
                alias_moved=False
                for moved_source in new_reaction_aliases:
                    if(original_alias in new_reaction_aliases[moved_source]):
                        alias_moved=True

                if(alias_moved==True):
                    if(original_source not in new_reaction_aliases):
                        new_reaction_aliases[original_source]=list()
                    new_reaction_aliases[original_source].append(original_alias)
                else:
                    if(original_source not in keep_reaction_aliases):
                        keep_reaction_aliases[original_source]=list()
                    if(original_alias not in keep_reaction_aliases[original_source]):
                        keep_reaction_aliases[original_source].append(original_alias)

    if(already_in_database is True):
#        print(original_rxn,keep_reaction_aliases)
#        print(disambig_rxn,new_reaction_aliases)
#        print("================================")
        pass

    reaction_aliases_dict[original_rxn]=keep_reaction_aliases
    reaction_aliases_dict[disambig_rxn]=new_reaction_aliases

    #Collect names and EC numbers for disambiguation
    ecname_disambig_dict={'from':{'id':original_rxn,
                                  'def':reactions_dict[original_rxn]['definition']},
                          'to':{'id':disambig_rxn,
                                'def':reactions_dict[disambig_rxn]['definition']},
                          'names':{},'ecs':{}}
    
    for name in reaction_names_dict[original_rxn]:
        ecname_disambig_dict['names'][name]="true"
    for ec in reaction_ecs_dict[original_rxn]:
        ecname_disambig_dict['ecs'][ec]="true"
    
    if(disambig_rxn in reaction_names_dict):
        for name in reaction_names_dict[original_rxn]:
            if(name in ecname_disambig_dict):
                ecname_disambig_dict['names'][name]="both"
            else:
                ecname_disambig_dict['names'][name]="false"
        for ec in reaction_ecs_dict[original_rxn]:
            if(ec in ecname_disambig_dict):
                ecname_disambig_dict['ecs'][ec]="both"
            else:
                ecname_disambig_dict['ecs'][ec]="false"

    disambiguation_object['reactions'].append(ecname_disambig_dict)

if(Dry_Run==0):
    print("Saving update to aliases")
    reactions_helper.saveAliases(reaction_aliases_dict)
    print("Saving disambiguated reactions")
    reactions_helper.saveReactions(reactions_dict)

if(len(disambiguation_object['reactions'])>0):
    print("Saving reaction names and ECs for disambiguation")
    with open('Objects/'+Old_Cpd_ID+'_'+metadata['user']+'_Reactions_Object.json','w') as f:
        json_string = json.dumps(disambiguation_object,
                                 indent=4, sort_keys=True,
                                 separators=(',', ': '), ensure_ascii=True)
        f.write(json_string)

#Scripts to run afterwards
#../Structures/List_ModelSEED_Structures.py
#../Structures/Update_Structure_Formula_Charge.py
#../Biochemistry/Update_Compound_Aliases.py
#../Biochemistry/Rebalance_Reactions.py (very important)
#../Biochemistry/Adjust_Reaction_Protons.py
#../Biochemistry/Adjust_Reaction_Water.py
#../Biochemistry/Merge_Reactions.py (merges may happen because of water)
