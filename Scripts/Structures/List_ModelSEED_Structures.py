#!/usr/bin/env python
import os, sys

sys.path.append('../../Libs/Python')
from BiochemPy import Compounds #, Reactions

#Load Compounds
CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()

#Load Curated Structures
Ignored_Structures=dict()
with open("../../Biochemistry/Structures/Curation/Ignore_Structures.txt") as ignore_file:
    for line in ignore_file.readlines():
        array=line.split('\t')
        Ignored_Structures[array[0]]=1
ignore_file.close()

#Load Structures and Aliases
Structures_Dict = CompoundsHelper.loadStructures(["SMILE","InChIKey","InChI"],["KEGG","MetaCyc"])
MS_Aliases_Dict =  CompoundsHelper.loadMSAliases(["KEGG","MetaCyc"])


unique_structs_file = open("../../Biochemistry/Structures/Unique_ModelSEED_Structures.txt",'w')
master_structs_file = open("../../Biochemistry/Structures/All_ModelSEED_Structures.txt",'w')
unique_structs_file.write("ID\tType\tAliases\tStructure\n")
for msid in sorted(MS_Aliases_Dict.keys()):

    #Build collection of all structures for the ModelSEED ID
    Structs = dict()
    for source in 'KEGG','MetaCyc':
        if(source not in MS_Aliases_Dict[msid].keys()):
            continue

        for struct_type in sorted(Structures_Dict.keys()):
            for external_id in sorted(MS_Aliases_Dict[msid][source]):
                if(external_id not in Structures_Dict[struct_type]):
                    continue

                for struct_stage in sorted(Structures_Dict[struct_type][external_id].keys()):
                    if(struct_type not in Structs):
                        Structs[struct_type]=dict()

                    if(struct_stage not in Structs[struct_type]):
                        Structs[struct_type][struct_stage]=dict()

                    for structure in sorted(Structures_Dict[struct_type][external_id][struct_stage].keys()):
                            
                        #Write to master
                        master_structs_file.write("\t".join([msid,struct_type,struct_stage,external_id,source,structure])+"\n")    

                        if(external_id in Ignored_Structures):
                            continue

                        if(structure not in Structs[struct_type][struct_stage]):
                            Structs[struct_type][struct_stage][structure]=dict()

                        Structs[struct_type][struct_stage][structure][external_id]=source
                        
    ms_structure = "null"
    ms_structure_type = "null"
    ms_external_ids = list()

    if("SMILE" in Structs):
        if("Charged" in Structs["SMILE"]):
            if(len(Structs["SMILE"]["Charged"])==1):
                structure = list(Structs["SMILE"]["Charged"].keys())[0]
                ms_structure = structure
                ms_structure_type = "SMILE"
                ms_external_ids = Structs["SMILE"]["Charged"][structure].keys()
            else:
                #Establish rules for checking/curating InChI strings
                pass

        elif("Original" in Structs["SMILE"]):
            if(len(Structs["SMILE"]["Original"])==1):
                structure = list(Structs["SMILE"]["Original"].keys())[0]
                ms_structure = structure
                ms_structure_type = "SMILE"
                ms_external_ids = Structs["SMILE"]["Original"][structure].keys()
            else:
                #Establish rules for checking/curating InChI strings
                pass
    else:
        pass

    if(ms_structure != "null"):
        unique_structs_file.write("\t".join([msid,ms_structure_type,";".join(sorted(ms_external_ids)),ms_structure])+"\n")
    else:
        if("KEGG" in MS_Aliases_Dict[msid]):
            #We have a problem where a structure is available for both KEGG and MetaCyc, but not only does conflict arise
            #The MetaCyc InChI over-rides the KEGG SMILE
            continue

    ms_structure="null"
    if("InChIKey" in Structs):

        #Default to structures charged by MarvinBeans
        if("Charged" in Structs["InChIKey"]):

            if(len(Structs["InChIKey"]["Charged"])==1):
                structure = list(Structs["InChIKey"]["Charged"].keys())[0]
                ms_structure = structure
                ms_structure_type = "InChIKey"
                ms_external_ids = Structs["InChIKey"]["Charged"][structure].keys()
            else:
                #Establish rules for checking/curating InChIKey strings
                pass

        elif("Original" in Structs["InChIKey"]):
            if(len(Structs["InChIKey"]["Original"])==1):
                structure = list(Structs["InChIKey"]["Original"].keys())[0]
                ms_structure = structure
                ms_structure_type = "InChIKey"
                ms_external_ids = Structs["InChIKey"]["Original"][structure].keys()
            else:
                #Establish rules for checking/curating InChIKey strings
                pass

    if(ms_structure != "null"):
        unique_structs_file.write("\t".join([msid,ms_structure_type,";".join(sorted(ms_external_ids)),ms_structure])+"\n")

    ms_structure="null"
    if("InChI" in Structs):

        #Default to structures charged by MarvinBeans
        if("Charged" in Structs["InChI"]):

            if(len(Structs["InChI"]["Charged"])==1):
                structure = list(Structs["InChI"]["Charged"].keys())[0]
                ms_structure = structure
                ms_structure_type = "InChI"
                ms_external_ids = Structs["InChI"]["Charged"][structure].keys()
            else:
                #Establish rules for checking/curating InChI strings
                pass

        elif("Original" in Structs["InChI"]):
            if(len(Structs["InChI"]["Original"])==1):
                structure = list(Structs["InChI"]["Original"].keys())[0]
                ms_structure = structure
                ms_structure_type = "InChI"
                ms_external_ids = Structs["InChI"]["Original"][structure].keys()
            else:
                #Establish rules for checking/curating InChI strings
                pass

    if(ms_structure != "null"):
        unique_structs_file.write("\t".join([msid,ms_structure_type,";".join(sorted(ms_external_ids)),ms_structure])+"\n")

unique_structs_file.close()
