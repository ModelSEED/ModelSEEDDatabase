#!/usr/bin/env python
import os, sys
from csv import DictReader
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()
Aliases_Dict = CompoundsHelper.loadMSAliases()
Names_Dict = CompoundsHelper.loadNames()

Source_Classes=dict()
reader = DictReader(open('../../Biochemistry/Aliases/Source_Classifiers.txt'), dialect='excel-tab')
for line in reader:
    if(line['Source Type'] not in Source_Classes):
        Source_Classes[line['Source Type']]=dict()
    Source_Classes[line['Source Type']][line['Source ID']]=1

for cpd in sorted(Compounds_Dict.keys()):
    if(cpd not in Aliases_Dict):
        continue

    Cpd_Aliases=dict()
    Alias_Count=0
    for source_type in 'Primary Database', 'Secondary Database', 'Published Model':
        for source in sorted(Aliases_Dict[cpd].keys()):
        
            if(len(Cpd_Aliases)>4):
                break

            if(source == "BiGG1"):
                continue

            if(source in Source_Classes[source_type] or source == "BiGG"):
                if(source not in Cpd_Aliases):
                    Cpd_Aliases[source]=dict()
                for alias in Aliases_Dict[cpd][source]:
                    Cpd_Aliases[source][alias]=1
                    Alias_Count+=1

    Alias_List=list()
    for source in sorted(Cpd_Aliases.keys()):
        source_line=source+":"+"|".join(sorted(Cpd_Aliases[source].keys()))
        Alias_List.append(source_line)
        
    if(cpd in Names_Dict):
        name_line="name:"+"|".join(sorted(Names_Dict[cpd]['name']))
        Alias_List.append(name_line)

    Alias_Line = ";".join(Alias_List)
    Compounds_Dict[cpd]['aliases']=Alias_Line

CompoundsHelper.saveCompounds(Compounds_Dict)
