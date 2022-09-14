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
        re.search('^rxn\d{5}$',arguments[0]) is None or \
        ( len(arguments) == 2 and re.search('^rxn\d{5}$',arguments[1]) is None)):
    print("Error: script must be initiated with one or two ModelSEED reaction identifier(s)")
    sys.exit()

disambiguating_rxn=arguments[0]
merging_rxn=None
if(len(arguments)==2):
    merging_rxn=arguments[1]
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
# Collect reaction data
#
##########################################################
sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

if(disambiguating_rxn not in reactions_dict):
    print("Error: reaction "+disambiguating_rxn+" is not found in the ModelSEED database")
    sys.exit()

if(reactions_dict[disambiguating_rxn]['is_obsolete']==1):
    print("Warning: reaction "+disambiguating_rxn+" is obsolete, consider using the non-obsolete version")

Disambiguation_Object['from']={'id':disambiguating_rxn,'aliases':{},'names':{},
                               'definition':reactions_dict[disambiguating_rxn]['definition']}

Aliases_Dict = reactions_helper.loadMSAliases()
Names_Dict = reactions_helper.loadNames()

#For reverse lookup
reverse_aliases_dict = dict()
for rxn in Aliases_Dict:
    for source in Aliases_Dict[rxn]:
        for alias in Aliases_Dict[rxn][source]:
            if(alias not in reverse_aliases_dict):
                reverse_aliases_dict[alias]=dict()
            if(source not in reverse_aliases_dict[alias]):
                reverse_aliases_dict[alias][source]=dict()
            reverse_aliases_dict[alias][source][rxn]=1

names_dict=dict()
for name in Names_Dict[disambiguating_rxn]:
    names_dict[name]=True
Disambiguation_Object['from']['names']=names_dict

aliases_dict=dict()
for source in Aliases_Dict[disambiguating_rxn]:
    for alias in Aliases_Dict[disambiguating_rxn][source]:
        if(alias not in aliases_dict):
            aliases_dict[alias]=dict()
        aliases_dict[alias][source]=True
Disambiguation_Object['from']['aliases']=aliases_dict

#It's possible that the structure and/or alias is associated with another reaction
#Need to find and report
other_reactions=list()
for rxn in Aliases_Dict.keys():
    if(rxn == disambiguating_rxn):
        continue

    for source in Aliases_Dict[rxn].keys():
        for alias in Aliases_Dict[rxn][source]:
            if(alias in aliases_dict and rxn not in other_reactions):
                other_reactions.append(rxn)

#Array of possible other reactions
Disambiguation_Object['to']=list()
if(merging_rxn is not None):
    Disambiguation_Object['to'].append({'id':merging_rxn,'name':reactions_dict[merging_rxn]['name'],
                                        'definition':reactions_dict[merging_rxn]['definition']})
elif(len(other_reactions)==0):
    Disambiguation_Object['to'].append({'id':None,'name':None,'definition':None})
else:
    for other_rxn in sorted(other_reactions):
        Disambiguation_Object['to'].append({'id':other_rxn,'name':reactions_dict[other_rxn]['name'],
                                            'definition':reactions_dict[other_rxn]['definition']})

##########################################################
#
# Print reaction object
#
##########################################################

filename=disambiguating_rxn
if(merging_rxn is not None):
    filename+='_'+merging_rxn
filename+='_'+user+'_Object.json'

with open('Objects/'+filename,'w') as f:
    json_string = json.dumps(Disambiguation_Object,
                             indent=4, sort_keys=True,
                             separators=(',', ': '), ensure_ascii=True)
    f.write(json_string)
