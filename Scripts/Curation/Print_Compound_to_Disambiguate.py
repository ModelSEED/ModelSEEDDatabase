#!/usr/bin/env python
import os, sys, subprocess, time, pprint, copy, re, json

Compound="cpd11665"
Disambiguation_Object = {'metadata':{},'from':{},'to':{}}

##########################################################
#
# Collect metadata
#
##########################################################

user=subprocess.check_output(['git','config','--global','--get','user.name'],universal_newlines=True).strip()
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

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()

Disambiguation_Object['from']={'id':Compound,'structures':{},'aliases':{},'names':{},
                               'formula':Compounds_Dict[Compound]['formula'],
                               'charge':Compounds_Dict[Compound]['charge']}

Aliases_Dict = CompoundsHelper.loadMSAliases()
Names_Dict = CompoundsHelper.loadNames()
Structures_Dict = CompoundsHelper.loadStructures(["InChI","SMILE"],["KEGG","MetaCyc"])

names_dict=dict()
for name in Names_Dict[Compound]:
    names_dict[name]=True
Disambiguation_Object['from']['names']=names_dict

aliases_dict=dict()
for source in Aliases_Dict[Compound]:
    for alias in Aliases_Dict[Compound][source]:
        if(alias not in aliases_dict):
            aliases_dict[alias]=dict()
        aliases_dict[alias][source]=True
Disambiguation_Object['from']['aliases']=aliases_dict

#It's possible that the alias is associated with another compound
#Need to find and report
Other_Compounds=list()
for cpd in Aliases_Dict.keys():
    if(cpd == Compound):
        continue

    for source in Aliases_Dict[cpd].keys():
        for alias in Aliases_Dict[cpd][source]:
            if(alias in aliases_dict and cpd not in Other_Compounds):
                Other_Compounds.append(cpd)

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

#Array of possible other compounds
Disambiguation_Object['to']=list()
if(len(Other_Compounds)==0):
    Disambiguation_Object['to'].append({'id':None,'name':None,'formula':None,'charge':None,'mass':None})
else:
    for cpd in sorted(Other_Compounds):
        Disambiguation_Object['to'].append({'id':cpd,'name':Compounds_Dict[cpd]['name'],
                                            'formula':Compounds_Dict[cpd]['formula'],
                                            'charge':Compounds_Dict[cpd]['charge'],
                                            'mass':Compounds_Dict[cpd]['mass']})

##########################################################
#
# Print compound object
#
##########################################################

with open('Objects/'+Compound+'_'+user+'_Object.json','w') as f:
    json_string = json.dumps(Disambiguation_Object,
                             indent=4, sort_keys=True,
                             separators=(',', ': '), ensure_ascii=True)
    f.write(json_string)
