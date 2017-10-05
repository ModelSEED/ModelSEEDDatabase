#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

CompoundsHelper = Compounds()

#Load Compounds
Compounds_Dict = CompoundsHelper.loadCompounds()

#Load Structures
Structures_Dict = CompoundsHelper.loadStructures()

#Load Aliases
MS_Aliases_Dict =  CompoundsHelper.loadMSAliases()

#Source_Aliases_Dict =  CompoundsHelper.loadSourceAliases()

for msid in MS_Aliases_Dict.keys():
    ms_structure = "null"

    #prioritizes KEGG structures (arbitrarily)
    for source in ["KEGG","MetaCyc"]:
        if(source not in MS_Aliases_Dict[msid]):
            continue

        #prioritizes source aliases in alphabetical order
        for external_id in sorted(MS_Aliases_Dict[msid][source]):
            if(ms_structure != "null"):
                break

            if(external_id in Structures_Dict['SMILE']):

                #prioritizes structures that originated from 'Charged' Mol file
                for struct_stage in ["Charged","Original"]:
                    if(ms_structure != "null"):
                        break

                    if(struct_stage not in Structures_Dict['SMILE'][external_id]):
                        continue

                    #prioritizes structures in alphabetical order
                    for structure in sorted(Structures_Dict['SMILE'][external_id][struct_stage]):
#                        print source,external_id,struct_stage
                        ms_structure=structure

    Compounds_Dict[msid]['smile']=ms_structure

for msid in Compounds_Dict.keys():
    if(Compounds_Dict[msid]['smile'] is None):
        Compounds_Dict[msid]['smile']="null"

CompoundsHelper.saveCompounds(Compounds_Dict)
