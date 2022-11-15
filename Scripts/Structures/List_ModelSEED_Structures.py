#!/usr/bin/env python
import os
import sys
import json
import glob

sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

#Load Compounds
CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()

Structures_Root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
Formulas_Dict=dict()
for source in "KEGG","MetaCyc": #,"ChEBI","Rhea":
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
for file in glob.glob(Structures_Root+"Curation/*.txt"):
    with open(file) as ignore_file:
        for line in ignore_file.readlines():
            array=line.split('\t')
            Ignored_Structures[array[0]]=1

#Load Structures and Aliases
Structures_Dict = CompoundsHelper.loadStructures(["SMILE","InChIKey","InChI"],["KEGG","MetaCyc"]) #,"ChEBI","Rhea"])
MS_Aliases_Dict =  CompoundsHelper.loadMSAliases(["KEGG","MetaCyc"]) #,"ChEBI","Rhea"])

master_structs_file = open(Structures_Root+"All_ModelSEED_Structures.txt",'w')
unique_structs_file = open(Structures_Root+"Unique_ModelSEED_Structures.txt",'w')
unique_structs_file.write("ID\tType\tAliases\tFormula\tCharge\tStructure\n")
structure_conflicts_file = open("Structure_Conflicts.txt",'w')
formula_conflicts_file = open("Formula_Conflicts.txt",'w')
for msid in sorted(MS_Aliases_Dict.keys()):

    #Build collection of all structures for the ModelSEED ID
    Structs = dict()
    Formulas=dict()
    for source in 'KEGG','MetaCyc': #,'ChEBI','Rhea':
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

        #If there are structural conflicts we will collect the ids and strings and establish rules
        #But if there isn't a structural conflict, we will use the one
        if(struct_conflict == 0):

            #In order to replicate, print SMILE, InChIKey, InChI in order
            for structure_type in "SMILE","InChIKey","InChI":
                if(structure_type not in Structs):
                    continue

                aliases=dict()
                
                # Even though there is one standardized InChI structure, there's a chance of multiple SMILES
                # So we "choose" a structure, but collect all aliases in case
                structure = sorted(Structs[structure_type][struct_stage].keys())[0]
                for struct in Structs[structure_type][struct_stage].keys():
                    for alias in Structs[structure_type][struct_stage][struct]:
                        aliases[alias]=1

                unique_structs_file.write("\t".join((msid,\
                                                         structure_type,\
                                                         ";".join(sorted(aliases)),\
                                                         formula_charge_dict['formula'],\
                                                         formula_charge_dict['charge'],\
                                                         structure))+"\n")

        else:
            struct_conflicts = dict()
            sources_structures=dict()
            #First decide on source and alias before printing, decision made with either InChI or SMILE alone
            for structure_type in "InChIKey","SMILE":
                if(structure_type not in Structs):
                    continue

                for structure in Structs[structure_type][struct_stage]:
                    for alias in Structs[structure_type][struct_stage][structure]:
                        source = Structs[structure_type][struct_stage][structure][alias]
                        if(structure not in struct_conflicts):
                            struct_conflicts[structure]=dict()
                        if(source not in struct_conflicts):
                            struct_conflicts[structure][source]=dict()
                        struct_conflicts[structure][source][alias]=1
                        if(source not in sources_structures):
                            sources_structures[source]=dict()
                        sources_structures[source][structure]=1

                if(len(struct_conflicts)>0):
                    break
            
            chosen_structure = None

            #Choose structure that is entirely identical from multiple sources
            #(indicates incorrect merger of additional id from same source)
            #At time of writing, there's two modelseed compounds that have more than
            #one structure with multiple sources, attempt to choose the ones that have
            #stereochemistry
            chosen_structures=dict()
            for structure in struct_conflicts:
                if(len(struct_conflicts[structure])>1):
                    chosen_structures[structure]=1

            if(len(chosen_structures)>0):
                if(len(chosen_structures)==1):
                    chosen_structure = list(chosen_structures.keys())[0]
                else:
                    for structure in chosen_structures:
                        #Avoid lack of stereochemistry, at time of writing, never happens
                        #For SMILE string
                        if('UHFFFAOYSA' not in structure):
                            chosen_structure = structure
                            break
                
            else:
                #Here, each different structure has a single source
                #Now, here, we can do no more for SMILE string, so pick MetaCyc if possible
                if(structure_type == "SMILE"):
                    if('MetaCyc' in sources_structures):
                        chosen_structure = sorted(sources_structures['MetaCyc'])[0]
                    elif('KEGG' in sources_structures):
                        chosen_structure = sorted(sources_structures['KEGG'])[0]
                    elif('ChEBI' in sources_structures):
                        chosen_structure = sorted(sources_structures['ChEBI'])[0]
                    elif('Rhea' in sources_structures):
                        chosen_structure = sorted(sources_structures['Rhea'])[0]

                else:
                    #So now we break down the InChIKey and find structures with the same connectivity
                    connected_structures = dict()
                    for structure in struct_conflicts:
                        connectivity = structure.split('-')[0]
                        if(connectivity not in connected_structures):
                            connected_structures[connectivity] = dict()
                        connected_structures[connectivity][structure]=1

                    chosen_connectivity = None
                    for connectivity in connected_structures:
                        if(len(connected_structures[connectivity])>1):
                            #At time of writing, only happens once per compound
                            chosen_connectivity = connectivity
                    
                    if(chosen_connectivity is not None):                        
                        #First pick structure that has stereo
                        stereo_structures = dict()
                        for structure in connected_structures[chosen_connectivity]:
                            if('UHFFFAOYSA' not in structure):
                                stereo_structures[structure]=1
                        if(len(stereo_structures)==1):
                            chosen_structure=list(stereo_structures.keys())[0]
                        else:
                            if('MetaCyc' in sources_structures):
                                chosen_structure = sorted(sources_structures['MetaCyc'])[0]
                            elif('KEGG' in sources_structures):
                                chosen_structure = sorted(sources_structures['KEGG'])[0]
                            elif('ChEBI' in sources_structures):
                                chosen_structure = sorted(sources_structures['ChEBI'])[0]
                            elif('Rhea' in sources_structures):
                                chosen_structure = sorted(sources_structures['Rhea'])[0]

                    if(chosen_structure is None):
                        #Here we have structures with the same formula, but different connectivity
                        #For now, we pick MetaCyc
                        if('MetaCyc' in sources_structures):
                            chosen_structure = sorted(sources_structures['MetaCyc'])[0]
                        elif('KEGG' in sources_structures):
                            chosen_structure = sorted(sources_structures['KEGG'])[0]
                        elif('ChEBI' in sources_structures):
                            chosen_structure = sorted(sources_structures['ChEBI'])[0]
                        elif('Rhea' in sources_structures):
                            chosen_structure = sorted(sources_structures['Rhea'])[0]

            #Now we have the chosen structure, we collect aliases of chosen structure
            #And use them to make sure we consistently use the right structure of each type
            chosen_aliases=dict()
            for source in struct_conflicts[chosen_structure]:
                for alias in struct_conflicts[chosen_structure][source]:
                    chosen_aliases[alias]=1

            structure_to_use = None
            for structure_type in "SMILE","InChIKey","InChI":
                if(structure_type not in Structs):
                    continue

                for alias in chosen_aliases:
                    for structure in Structs[structure_type][struct_stage]:
                        if(alias in Structs[structure_type][struct_stage][structure]):
                            structure_to_use=structure

                #Now that we've determined the structure to use for that type
                #We collect all aliases regardless of structure
                aliases=dict()
                for structure in Structs[structure_type][struct_stage]:
                    for alias in Structs[structure_type][struct_stage][structure]:
                        aliases[alias]=1

                #Finally, write to file
                unique_structs_file.write("\t".join((msid,\
                                                         structure_type,\
                                                         ";".join(sorted(aliases)),\
                                                         formula_charge_dict['formula'],\
                                                         formula_charge_dict['charge'],\
                                                         structure_to_use))+"\n")

    if(struct_conflict==1):
        for structure in Structs[struct_type][struct_stage]:
            for external_id in Structs[struct_type][struct_stage][structure]:
                structure_conflicts_file.write("\t".join((msid,struct_type,struct_stage,structure,external_id,
                                                          Structs[struct_type][struct_stage][structure][external_id]))+"\n")

    if(formula_conflict==1):
        for formula in Formulas[struct_type][struct_stage]:
            for external_id in Formulas[struct_type][struct_stage][formula]:
                formula_dict=json.loads(formula)
                formula_conflicts_file.write("\t".join((msid,struct_type,struct_stage,formula_dict['formula'],formula_dict['charge'],external_id,
                                                          Formulas[struct_type][struct_stage][formula][external_id]))+"\n")
