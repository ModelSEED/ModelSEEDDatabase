#!/usr/bin/env python
import os
import sys
import subprocess
import time
import copy
import re
import json

arguments = list(sys.argv)
#Pop filename
arguments = arguments[1:]
if( len(arguments) < 1 or len(arguments) > 2 or \
        re.search('^cpd\d{5}$',arguments[0]) is None or \
        ( len(arguments) == 2 and re.search('^cpd\d{5}$',arguments[1]) is None)):
    print("Error: script must be initiated with one or two ModelSEED compound identifier(s)")
    sys.exit()

disambiguating_cpd=arguments[0]
merging_cpd=None
if(len(arguments)==2):
    merging_cpd=arguments[1]
Disambiguation_Object = {'metadata':{},'from':{},'to':{},'method':"disambiguation"}

##########################################################
#
# Collect metadata
#
##########################################################

user=subprocess.check_output(['git','config','--global','--get','user.name'],universal_newlines=True).strip()
user=re.sub(' ','_',user)
Disambiguation_Object['metadata']['user']=user

output=subprocess.check_output(['git','remote','show','origin'],universal_newlines=True)
remote=""
for line in output.split('\n'):
    if('Fetch URL' in line):
        remote=line.split(' ')[-1]
Disambiguation_Object['metadata']['remote']=remote

output=subprocess.check_output(['git','rev-parse','HEAD'],universal_newlines=True)
commit=output.strip()
Disambiguation_Object['metadata']['commit']=commit

output=subprocess.check_output(['git','rev-parse','--abbrev-ref','HEAD'],universal_newlines=True)
branch=output.strip()
Disambiguation_Object['metadata']['branch']=branch

time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
Disambiguation_Object['metadata']['date_time']=time_str

##########################################################
#
# Collect compound data
#
##########################################################
from BiochemPy import Reactions, Compounds

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

if(disambiguating_cpd not in compounds_dict):
    print("Error: compound "+disambiguating_cpd+" is not found in the ModelSEED database")
    sys.exit()

if(compounds_dict[disambiguating_cpd]['is_obsolete']==1):
    print("Warning: compound "+disambiguating_cpd+" is obsolete, consider using the non-obsolete version")

Disambiguation_Object['from']={'id':disambiguating_cpd,'structures':{},'aliases':{},'names':{},
                               'formula':compounds_dict[disambiguating_cpd]['formula'],
                               'charge':compounds_dict[disambiguating_cpd]['charge'],
                               'mass':compounds_dict[disambiguating_cpd]['mass']}

Aliases_Dict = compounds_helper.loadMSAliases()
Names_Dict = compounds_helper.loadNames()
Structures_Dict = compounds_helper.loadStructures(["InChI","SMILE"],["KEGG","MetaCyc"])

#For reverse lookup
reverse_aliases_dict = dict()
for cpd in Aliases_Dict:
    for source in Aliases_Dict[cpd]:
        for alias in Aliases_Dict[cpd][source]:
            if(alias not in reverse_aliases_dict):
                reverse_aliases_dict[alias]=dict()
            if(source not in reverse_aliases_dict[alias]):
                reverse_aliases_dict[alias][source]=dict()
            reverse_aliases_dict[alias][source][cpd]=1

reverse_structures_dict = dict()
for type in Structures_Dict:
    for alias in Structures_Dict[type]:
        for stage in Structures_Dict[type][alias]:
            for structure in Structures_Dict[type][alias][stage]:
                if(structure not in reverse_structures_dict):
                    reverse_structures_dict[structure]=dict()
                if(stage not in reverse_structures_dict[structure]):
                    reverse_structures_dict[structure][stage]=dict()
                reverse_structures_dict[structure][stage][alias]=1

names_dict=dict()
for name in Names_Dict[disambiguating_cpd]:
    names_dict[name]=True
Disambiguation_Object['from']['names']=names_dict

aliases_dict=dict()
for source in Aliases_Dict[disambiguating_cpd]:
    for alias in Aliases_Dict[disambiguating_cpd][source]:
        if(alias not in aliases_dict):
            aliases_dict[alias]=dict()
        aliases_dict[alias][source]=True
Disambiguation_Object['from']['aliases']=aliases_dict

#It's possible that the structure and/or alias is associated with another compound
#Need to find and report
other_compounds=list()
for cpd in Aliases_Dict.keys():
    if(cpd == disambiguating_cpd):
        continue

    for source in Aliases_Dict[cpd].keys():
        for alias in Aliases_Dict[cpd][source]:
            if(alias in aliases_dict and cpd not in other_compounds):
                other_compounds.append(cpd)

structures_dict=dict()
inchi=False
for alias in Structures_Dict['InChI'].keys():
    if(alias in aliases_dict):
        inchi=True
        if(alias not in structures_dict):
            structures_dict[alias]=dict()
        if('Charged' in Structures_Dict['InChI'][alias]):
            for structure in Structures_Dict['InChI'][alias]['Charged']:
                structures_dict[alias][structure]=True
        elif('Original' in Structures_Dict['InChI'][alias]):
            for structure in Structures_Dict['InChI'][alias]['Original']:
                structures_dict[alias][structure]=True

if(inchi is not True):
    for alias in Structures_Dict['SMILE'].keys():
        if(alias in aliases_dict):
            if(alias not in structures_dict):
                structures_dict[alias]=dict()
            if('Charged' in Structures_Dict['SMILE'][alias]):
                for structure in Structures_Dict['SMILE'][alias]['Charged']:
                    structures_dict[alias][structure]=True
            elif('Original' in Structures_Dict['SMILE'][alias]):
                for structure in Structures_Dict['SMILE'][alias]['Original']:
                    structures_dict[alias][structure]=True
Disambiguation_Object['from']['structures']=structures_dict

#Find other_compounds using structures
#This should work regardless of whether the structure is an InChI or a SMILE
for current_alias in structures_dict:
    for current_structure in structures_dict[current_alias]:
        if(current_structure in reverse_structures_dict):
            if('Charged' in reverse_structures_dict[current_structure]):
                for other_alias in reverse_structures_dict[current_structure]['Charged']:
                    for source in reverse_aliases_dict[other_alias]:
                        for cpd in reverse_aliases_dict[other_alias][source]:
                            if(cpd != disambiguating_cpd and cpd not in other_compounds):
                                other_compounds.append(cpd)

#Array of possible other compounds
Disambiguation_Object['to']=list()
#Here use specific compound listed by user
if(merging_cpd is not None):
    Disambiguation_Object['to'].append({'id':merging_cpd,'name':compounds_dict[merging_cpd]['name'],
                                        'formula':compounds_dict[merging_cpd]['formula'],
                                        'charge':compounds_dict[merging_cpd]['charge'],
                                        'mass':compounds_dict[merging_cpd]['mass']})
elif(len(other_compounds)==0):
    Disambiguation_Object['to'].append({'id':None,'name':None,'formula':None,'charge':"0",'mass':"0"})
else:
    for other_cpd in sorted(other_compounds):
        Disambiguation_Object['to'].append({'id':other_cpd,'name':compounds_dict[other_cpd]['name'],
                                            'formula':compounds_dict[other_cpd]['formula'],
                                            'charge':compounds_dict[other_cpd]['charge'],
                                            'mass':compounds_dict[other_cpd]['mass']})

##########################################################
#
# Print compound object
#
##########################################################

filename=disambiguating_cpd
if(merging_cpd is not None):
    filename+='_'+merging_cpd
filename+='_'+user+'_Object.json'

with open('Objects/'+filename,'w') as f:
    json_string = json.dumps(Disambiguation_Object,
                             indent=4, sort_keys=True,
                             separators=(',', ': '), ensure_ascii=True)
    f.write(json_string)
