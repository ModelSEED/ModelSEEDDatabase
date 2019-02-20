#!/usr/bin/env python
import os
import sys
import json
from BiochemPy import Compounds

#Load Compounds
CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()

Structures_Root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
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
                    Formulas_Dict[source][struct_type][struct_stage][array[0]]={'formula':array[1],'charge':array[2]}

#Load Curated Structures
Ignored_Structures=dict()
#with open(Structures_Root+"Curation/Ignore_Structures.txt") as ignore_file:
#    for line in ignore_file.readlines():
#        array=line.split('\t')
#        Ignored_Structures[array[0]]=1
#ignore_file.close()

#Load Structures and Aliases
Structures_Dict = CompoundsHelper.loadStructures(["SMILE","InChIKey","InChI"],["KEGG","MetaCyc"])
MS_Aliases_Dict =  CompoundsHelper.loadMSAliases(["KEGG","MetaCyc"])

master_structs_file = open(Structures_Root+"All_ModelSEED_Structures.txt",'w')
unique_structs_file = open(Structures_Root+"Unique_ModelSEED_Structures.txt",'w')
unique_structs_file.write("ID\tType\tAliases\tFormula\tCharge\tStructure\n")
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

                        formula_charge_dict={'formula':"null",'charge':"null"}

                        if(struct_type in Formulas_Dict[source] and external_id in Formulas_Dict[source][struct_type][struct_stage]):
                            formula_charge_dict = Formulas_Dict[source][struct_type][struct_stage][external_id]

                        #Write to master
                        master_structs_file.write("\t".join([msid,struct_type,struct_stage,external_id,source,\
                                                                 formula_charge_dict['formula'],\
                                                                 formula_charge_dict['charge'],\
                                                                 structure])+"\n")

                        if(external_id in Ignored_Structures):
                            continue

                        if(structure not in Structs[struct_type][struct_stage]):
                            Structs[struct_type][struct_stage][structure]=dict()
                        Structs[struct_type][struct_stage][structure][external_id]=source

                        formula_charge_json = json.dumps(formula_charge_dict)
                        if(formula_charge_json not in Formulas[struct_type][struct_stage]):
                            Formulas[struct_type][struct_stage][formula_charge_json]=dict()
                        Formulas[struct_type][struct_stage][formula_charge_json][external_id]=source

    if(len(Structs.keys())==0):
        continue

    #Priority is:
    #Charged InChI
    #Original InChI
    #Charged SMILE
    #Original SMILE

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
        print("Warning: no structures used for "+msid)
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
        #Only one formula/charge combination possible here
        formula_charge_dict=json.loads(list(Formulas[struct_type][struct_stage].keys())[0])

        #In order to replicate, print SMILE, InChIKey, InChI in order
        for structure_type in "SMILE","InChIKey","InChI":
            if(structure_type not in Structs):
                continue

            #It's possible that there's a structural conflict but formula/charge combination is unique
            #In which case, we'd still like to use a single structure, and we will default to MetaCyc
            aliases=dict()
            sources_structures=dict()
            structure = None
            for structure in sorted(Structs[structure_type][struct_stage].keys()):
                for alias in Structs[structure_type][struct_stage][structure]:
                    aliases[alias]=1
                    source = Structs[structure_type][struct_stage][structure][alias]
                    if(source not in sources_structures):
                        sources_structures[source]=dict()
                    sources_structures[source][structure]=1

            if("MetaCyc" in sources_structures):
                structure = sorted(list(sources_structures["MetaCyc"].keys()))[0]
            elif("KEGG" in sources_structures):
                structure = sorted(list(sources_structures["KEGG"].keys()))[0]

            if(structure is None):
                #At time of writing this never happens
                continue

            unique_structs_file.write("\t".join((msid,\
                                                     structure_type,\
                                                     ";".join(sorted(aliases)),\
                                                     formula_charge_dict['formula'],\
                                                     formula_charge_dict['charge'],\
                                                     structure))+"\n")

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
