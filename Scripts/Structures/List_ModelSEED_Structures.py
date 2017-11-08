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
Structures_Dict = CompoundsHelper.loadStructures(["SMILE","InChIKey","InChI"],["KEGG","MetaCyc"])

#Load Aliases
MS_Aliases_Dict =  CompoundsHelper.loadMSAliases(["KEGG","MetaCyc"])

#Source_Aliases_Dict =  CompoundsHelper.loadSourceAliases()

ms_structs_file = open("../../Biochemistry/Structures/ModelSEED/ModelSEED_Strings.txt",'w')
for msid in sorted(MS_Aliases_Dict.keys()):
    #Build collection of all structures for the ModelSEED ID
    Structs = dict()
    for source in MS_Aliases_Dict[msid]:
        for struct_type in Structures_Dict.keys():
            for external_id in sorted(MS_Aliases_Dict[msid][source]):
                if(external_id not in Structures_Dict[struct_type]):
                    continue

                for struct_stage in Structures_Dict[struct_type][external_id].keys():
                    if(struct_type not in Structs):
                        Structs[struct_type]=dict()


                    if(struct_stage not in Structs[struct_type]):
                        Structs[struct_type][struct_stage]=dict()

                    for structure in Structures_Dict[struct_type][external_id][struct_stage].keys():
                        if(structure not in Structs[struct_type][struct_stage]):
                            Structs[struct_type][struct_stage][structure]=dict()
                            
                        Structs[struct_type][struct_stage][structure][external_id]=source

    ms_structure = "null"
    ms_structure_type = "null"
    ms_external_ids = list()

    if("InChI" in Structs):

        #Default to structures charged by MarvinBeans
        if("Charged" in Structs["InChI"]):

            if(len(Structs["InChI"]["Charged"])==1):
                structure = Structs["InChI"]["Charged"].keys()[0]
                ms_structure = structure
                ms_structure_type = "InChI"
                ms_external_ids = Structs["InChI"]["Charged"][structure].keys()
            else:
                #Establish rules for checking/curating InChI strings
                pass

        elif("Original" in Structs["InChI"]):
            if(len(Structs["InChI"]["Original"])==1):
                structure = Structs["InChI"]["Original"].keys()[0]
                ms_structure = structure
                ms_structure_type = "InChI"
                ms_external_ids = Structs["InChI"]["Original"][structure].keys()
            else:
                #Establish rules for checking/curating InChI strings
                pass
 
    elif("SMILE" in Structs):
        if(len(Structs["SMILE"]["Original"])==1):
            structure = Structs["SMILE"]["Original"].keys()[0]
            ms_structure = structure
            ms_structure_type = "SMILE"
            ms_external_ids = Structs["SMILE"]["Original"][structure].keys()
        else:
            #Establish rules for checking/curating InChI strings
            pass
    else:
        pass

    if(ms_structure == "null"):
        continue

    ms_structs_file.write("\t".join([msid,ms_structure_type,";".join(sorted(ms_external_ids)),ms_structure])+"\n")

ms_structs_file.close()
