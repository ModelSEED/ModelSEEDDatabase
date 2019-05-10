#!/usr/bin/env python
import re
from BiochemPy import Reactions
from csv import DictReader

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()
Aliases_Dict = ReactionsHelper.loadMSAliases()
Names_Dict = ReactionsHelper.loadNames()
ECs_Dict = ReactionsHelper.loadECs()
Pwys_Dict = ReactionsHelper.loadPathways()

Source_Classes=dict()
reader = DictReader(open('../../Biochemistry/Aliases/Source_Classifiers.txt'), dialect='excel-tab')
for line in reader:
    if(line['Source Type'] not in Source_Classes):
        Source_Classes[line['Source Type']]=dict()
    Source_Classes[line['Source Type']][line['Source ID']]=1

for rxn in sorted(Reactions_Dict.keys()):
    if(rxn not in Aliases_Dict):
        continue

    Rxn_Source_Aliases=dict()
    Rxn_Aliases=dict()
    Alias_Count=0
    for source_type in 'Primary Database', 'Secondary Database', 'Published Model':
        for source in sorted(Aliases_Dict[rxn].keys()):
        
            if(len(Rxn_Source_Aliases)>4):
                break

            if(source == "BiGG1"):
                continue

            if("KEGG" in source and len(source)>4):
                continue

            if(source in Source_Classes[source_type] or source == "BiGG"):
                if(source not in Rxn_Source_Aliases):
                    Rxn_Source_Aliases[source]=dict()
                for alias in Aliases_Dict[rxn][source]:
                    Rxn_Aliases[alias]=1

                    if("Cyc" in source and bool(re.search('.[a-z]$',alias))):
                        alias=re.sub(".[a-z]$","",alias)
                    Rxn_Source_Aliases[source][alias]=1
                    Rxn_Aliases[alias]=1

    Alias_List=list()
    for source in sorted(Rxn_Source_Aliases.keys()):
        source_line=source+": "+"; ".join(sorted(Rxn_Source_Aliases[source].keys()))
        Alias_List.append(source_line)
        
    if(rxn in Names_Dict):
        name_list=list()
        for name in Names_Dict[rxn]:
            if(name in Rxn_Aliases):
                continue

            if(bool(re.search('.[a-z]$',name)) and re.sub(".[a-z]$","",name) in Rxn_Aliases):
                continue

            name_list.append(name)
            
        name_line="Name: "+"; ".join(sorted(name_list))
        Alias_List.append(name_line)

    if(len(Alias_List)==0):
        Reactions_Dict[rxn]['aliases']="null"
    else:
        Reactions_Dict[rxn]['aliases']=Alias_List

    ECs_List=list()
    if(rxn in ECs_Dict):
        for ec in sorted(ECs_Dict[rxn]):
            Rxn_Aliases[ec]=1
            ECs_List.append(ec)

    if(len(ECs_List)==0):
        Reactions_Dict[rxn]['ec_numbers']="null"
    else:
        Reactions_Dict[rxn]['ec_numbers']=ECs_List

    Pwys_List=list()
    if(rxn in Pwys_Dict):
        for biochem in Pwys_Dict[rxn]:
            pwy_list=list()
            for pwy in Pwys_Dict[rxn][biochem]:
                pwy_list.append(pwy)
            pwy_string = biochem+": "+"; ".join(pwy_list)
            Pwys_List.append(pwy_string)

    if(len(Pwys_List)==0):
        Reactions_Dict[rxn]['pathways']="null"
    else:
        Reactions_Dict[rxn]['pathways']=Pwys_List

ReactionsHelper.saveReactions(Reactions_Dict)
