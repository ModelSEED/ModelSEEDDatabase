#!/usr/bin/env python
import os, sys
from difflib import Differ

sys.path.append('../../Libs/Python')
from BiochemPy import Compounds #, Reactions

DifferObj=Differ()

#Load Compounds
CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()

Structures_Root="../../Biochemistry/Structures/"
Formulas_Dict=dict()
for source in "KEGG","MetaCyc":
    if(source not in Formulas_Dict):
        Formulas_Dict[source]=dict()

    for struct_type in "InChI","SMILE":
        if(struct_type not in Formulas_Dict[source]):
            Formulas_Dict[source][struct_type]=dict()

        for struct_stage in "Charged","Original":
            if(struct_stage not in Formulas_Dict[source][struct_type]):
                Formulas_Dict[source][struct_type][struct_stage]=dict()

            file_name=Structures_Root+source+"/"+struct_type+"_"+struct_stage+"_Formulas_Charges.txt"
            with open(file_name) as file_handle:
                for line in file_handle.readlines():
                    line=line.strip()
                    array=line.split('\t')
                    Formulas_Dict[source][struct_type][struct_stage][array[0]]=array[1]+array[2]

#Load Curated Structures
Ignored_Structures=dict()
with open(Structures_Root+"Curation/Ignore_Structures.txt") as ignore_file:
    for line in ignore_file.readlines():
        array=line.split('\t')
#        Ignored_Structures[array[0]]=1
ignore_file.close()

#Load Structures and Aliases
Structures_Dict = CompoundsHelper.loadStructures(["SMILE","InChIKey","InChI"],["KEGG","MetaCyc"])
MS_Aliases_Dict =  CompoundsHelper.loadMSAliases(["KEGG","MetaCyc"])

master_structs_file = open("../../Biochemistry/Structures/All_ModelSEED_Structures.txt",'w')
unique_structs_file = open("../../Biochemistry/Structures/Unique_ModelSEED_Structures.txt",'w')
unique_structs_file.write("ID\tType\tAliases\tStructure\n")
structure_conflicts_file = open("Structure_Conflicts.txt",'w')
formula_conflicts_file = open("Formula_Conflicts.txt",'w')
for msid in sorted(MS_Aliases_Dict.keys()):

    #Build collection of all structures for the ModelSEED ID
    Structs = dict()
    Formulas=dict()
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
                        Formulas[struct_type]=dict()

                    if(struct_stage not in Structs[struct_type]):
                        Structs[struct_type][struct_stage]=dict()
                        Formulas[struct_type][struct_stage]=dict()

                    for structure in sorted(Structures_Dict[struct_type][external_id][struct_stage].keys()):
                            
                        #Write to master
                        master_structs_file.write("\t".join([msid,struct_type,struct_stage,external_id,source,structure])+"\n")    

                        if(external_id in Ignored_Structures):
                            continue

                        if(structure not in Structs[struct_type][struct_stage]):
                            Structs[struct_type][struct_stage][structure]=dict()
                        Structs[struct_type][struct_stage][structure][external_id]=source

                        if(struct_type in Formulas_Dict[source]):

                            #There's a chance that extracting formula failed, see CHLOROPHYLLIDE-A
                            if(external_id not in Formulas_Dict[source][struct_type][struct_stage]):
                                continue

                            formula = Formulas_Dict[source][struct_type][struct_stage][external_id]
                            if(formula not in Formulas[struct_type][struct_stage]):
                                Formulas[struct_type][struct_stage][formula]=dict()
                            Formulas[struct_type][struct_stage][formula][external_id]=source

    if(len(Structs.keys())==0):
        continue

    #Charged InChI
    #Single Formula? Remember!
    #Check to see if, for one formula, there are multiple charges
    #If single formula, then print all InChI, InChIKey, SMILE (select SMILE with same formula? select SMILE from KEGG?)
    #If not single formula, print conflict
    #If not Charged InChI
    #Check Original InChI (repeat)
    #Check Original SMILE (repeat)

    struct_type=None
    struct_stage=None
    if("InChI" in Structs):
        struct_type="InChI"
        if("Charged" in Structs[struct_type]):
            struct_stage="Charged"
        elif("Original" in Structs[struct_type]):
            struct_stage="Original"
    elif("SMILE" in Structs):
        struct_type="SMILE"
        if("Charged" in Structs[struct_type]):
            struct_stage="Charged"
        elif("Original" in Structs[struct_type]):
            struct_stage="Original"

    if(struct_type is None or struct_stage is None):
        #At time of writing, this doesn't happen
        continue

    struct_pass=0
    struct_conflict=0
    formula_conflict=0
    if(len(Structs[struct_type][struct_stage].keys())==1):
        struct_pass=1
    elif(len(Formulas[struct_type][struct_stage].keys())==1):
        struct_conflict=1
        struct_pass=1
    else:
        struct_conflict=1
        formula_conflict=1
        pass

    if(struct_pass):
        #In order to replicate, print SMILE, InChIKey, InChI in order
        for structure_type in "SMILE","InChIKey","InChI":
            if(structure_type not in Structs):
                continue
            for structure in sorted(Structs[structure_type][struct_stage].keys()):
                unique_structs_file.write("\t".join((msid,structure_type,";".join(sorted(Structs[structure_type][struct_stage][structure])),structure))+"\n")

    if(struct_conflict==1):
        for structure in Structs[struct_type][struct_stage]:
            for external_id in Structs[struct_type][struct_stage][structure]:
                structure_conflicts_file.write("\t".join((msid,struct_type,struct_stage,structure,external_id,
                                                          Structs[struct_type][struct_stage][structure][external_id]))+"\n")

    if(formula_conflict==1):
        for formula in Formulas[struct_type][struct_stage]:
            for external_id in Formulas[struct_type][struct_stage][formula]:
                formula_conflicts_file.write("\t".join((msid,struct_type,struct_stage,formula,external_id,
                                                          Formulas[struct_type][struct_stage][formula][external_id]))+"\n")
