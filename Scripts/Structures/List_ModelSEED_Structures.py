#!/usr/bin/env python
import os
import sys
import json
import glob

#################################################################
## Load Compound Objects into memory
#################################################################
sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

#Load Compounds
CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()

#################################################################
## Load Formula Strings from file
#################################################################

Structures_Root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
Formulas_Dict=dict()
for source in "KEGG","MetaCyc","ChEBI","Rhea":
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

#################################################################
## Load Curated Picks for Structures
#################################################################

#Load Curated Structures
Ignored_Structures=dict()
for file in glob.glob(Structures_Root+"Curation/*.txt"):
    with open(file) as ignore_file:
        for line in ignore_file.readlines():
            array=line.split('\t')
            Ignored_Structures[array[0]]=1

#################################################################
## Load Structures and Aliases
#################################################################

#Load Structures and Aliases
Structures_Dict = CompoundsHelper.loadStructures(["SMILE","InChIKey","InChI"],["KEGG","MetaCyc","ChEBI","Rhea"])
MS_Aliases_Dict =  CompoundsHelper.loadMSAliases(["KEGG","MetaCyc","ChEBI","Rhea"])

#################################################################
## Open filehandles for writing
#################################################################

master_structs_file = open(Structures_Root+"All_ModelSEED_Structures.txt",'w')
unique_structs_file = open(Structures_Root+"Unique_ModelSEED_Structures.txt",'w')
unique_structs_file.write("ID\tType\tAliases\tFormula\tCharge\tStructure\n")
structure_conflicts_file = open("Structure_Conflicts.txt",'w')
formula_conflicts_file = open("Formula_Conflicts.txt",'w')

#################################################################
## Iterate through ModelSEED identifiers
#################################################################

for msid in sorted(MS_Aliases_Dict.keys()):

    #################################################################
    ## For the ModelSEED compound build dict of all structures
    #################################################################

    Structs = dict()
    Formulas=dict()
    for source in 'KEGG','MetaCyc','ChEBI','Rhea':
        if(source not in MS_Aliases_Dict[msid].keys()):
            continue

        #################################################################
        ## Iterate through types, sources, ids
        #################################################################

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

                        #################################################################
                        ## Write all structures to master 'All_ModelSEED_Strutures.txt'
                        #################################################################

                        master_structs_file.write("\t".join([msid,struct_type,struct_stage,external_id,source,\
                                                                 formula_charge_dict['formula'],\
                                                                 formula_charge_dict['charge'],\
                                                                 structure])+"\n")
                        
                        #################################################################
                        ## Skip if curated structure designated to be ignored
                        #################################################################

                        if(external_id in Ignored_Structures):
                            continue

                        #################################################################
                        ## Populate structures dictionary
                        #################################################################

                        if(structure not in Structs[struct_type][struct_stage]):
                            Structs[struct_type][struct_stage][structure]=dict()
                        Structs[struct_type][struct_stage][structure][external_id]=source

                        formula_charge_json = json.dumps(formula_charge_dict)
                        if(formula_charge_json not in Formulas[struct_type][struct_stage]):
                            Formulas[struct_type][struct_stage][formula_charge_json]=dict()
                        Formulas[struct_type][struct_stage][formula_charge_json][external_id]=source
                        
    #################################################################
    ## Skip if no structures were collected
    #################################################################

    if(len(Structs.keys())==0):
        continue
    
    #################################################################
    ## Prioritized which type and stage for the structure for comparison
    ## Priority Order is:
    ## 1) Charged InChI
    ## 2) Original InChI
    ## 3) Charged SMILE
    ## 4) Original SMILE
    #################################################################

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
                            
    #################################################################
    ## Now the prioritized type and stage has been established
    ## We look to see whether or not they have the same structure
    ## from different sources
    ##
    ## If there's only one structure string, or there's multiple
    ## structure strings that have the same formula, the structure
    ## Is considered a 'pass' NB: This will be re-written
    ## To accomodate for differing protonation states being depicted
    ## as charges
    #################################################################

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
                            
    #################################################################
    ## Here on, we consider the structure(s) in general for a
    ## ModelSEED compound to be unique to that compound though
    ## they may vary if there's still a conflict.
    ## The code in this section describes how we
    ## might 'pick' a structure.
    #################################################################

    if(struct_pass):

        #################################################################
        ## The formula is considered to be identical across structures
        ## But we've recently found this isn't always the case with some
        ## protonated alcohol groups
        #################################################################

        #Only one formula/charge combination possible here
        formula_charge_dict=json.loads(list(Formulas[struct_type][struct_stage].keys())[0])
        
        #################################################################
        ## If there are no structural conflicts then all the structures
        ## are identical and we pick one
        #################################################################

        if(struct_conflict == 0):

            #################################################################
            ## Iterate through the types for writing to file
            ## Keeping order consistent but it doesn't matter which
            #################################################################

            for structure_type in "SMILE","InChIKey","InChI":
                if(structure_type not in Structs):
                    continue

                #################################################################
                ## Even though there is one standardized InChI structure
                ## There can be more than one SMILE so here we choose one 
                #################################################################

                structure = sorted(Structs[structure_type][struct_stage].keys())[0]

                #################################################################
                ## But we collect all aliases for that type and stage
                ## This means that the aliases for different SMILE strings
                ## will be grouped together under a single SMILE string
                #################################################################
                
                aliases=dict()
                for struct in Structs[structure_type][struct_stage].keys():
                    for alias in Structs[structure_type][struct_stage][struct]:
                        aliases[alias]=1

                #################################################################
                ## We write them to file
                #################################################################

                unique_structs_file.write("\t".join((msid,\
                                                     structure_type,\
                                                     ";".join(sorted(aliases)),\
                                                     formula_charge_dict['formula'],\
                                                     formula_charge_dict['charge'],\
                                                     structure))+"\n")

        else:

            #################################################################
            ## Now here we have similar structures that are identical in
            ## formula but differ in connectivity somehow. The code in
            ## this section is about 'selecting' the primary structure from
            ## a database and an alias that represents the compound
            #################################################################
            
            struct_conflicts = dict()
            sources_structures=dict()

            #################################################################
            ## First we pick InChIKey over SMILE as the primary structure, if
            ## we can. Then we collect that type of structure
            #################################################################

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


                #################################################################
                ## Break if collected InChIKey structures
                #################################################################
                        
                if(len(struct_conflicts)>0):
                    break

            #################################################################
            ## Here we go through the structures and determine if any come
            ## from more than database to prioritize, and then if not, we
            ## break them down to see how their connectivity matches and then
            ## we prioritize certain databases.
            #################################################################

            chosen_structure = None

            #################################################################
            ## Find structures that are identical in more than one database
            #################################################################

            chosen_structures = dict()
            for structure in struct_conflicts:
                if(len(struct_conflicts[structure])>1):
                    chosen_structures[structure]=1

            if(len(chosen_structures)>0):
                if(len(chosen_structures)==1):

                        #################################################################
                        ## If there is one structure that is identical in more than one
                        ## database we pick that structure to be the primary one
                        #################################################################

                        chosen_structure = list(chosen_structures.keys())[0]
                    
                else:
                    for structure in chosen_structures:
                        
                        #################################################################
                        ## If more than one identical structures are found in more than
                        ## one database then we have to arbitrarily pick one
                        ## There's not many at all
                        #################################################################

                        #Avoid lack of stereochemistry, at time of writing, never happens
                        #For SMILE string
                        if('UHFFFAOYSA' not in structure):
                            chosen_structure = structure
                            break

                    if(chosen_structure is None):
                        
                        print(msid,chosen_structures)

                        chosen_structure = list(chosen_structures.keys())[0]
                        
            else:

                #################################################################
                ## Here, each different structure comes from one database so
                ## we have to pick one arbitrarily
                #################################################################

                #################################################################
                ## First, if it's a SMILE string and not InChI, it's difficult
                ## to parse, so just prioritize which database it came from
                ## and pick one
                #################################################################

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

                    #################################################################
                    ## Secondly, if it's not a SMILE, then as InChI, we can break
                    ## the structure down and compare its connectivity which is
                    ## the first InChIKey 'layer'. But, we need to re-configure
                    ## this to use InChI
                    #################################################################

                    connected_structures = dict()
                    for structure in struct_conflicts:
                        connectivity = structure.split('-')[0]
                        if(connectivity not in connected_structures):
                            connected_structures[connectivity] = dict()
                        connected_structures[connectivity][structure]=1
                        
                    #################################################################
                    ## Having split the InChiKey into its connectivity layers and
                    ## grouping them, we take the connectivity that has more than
                    ## one structure (i.e. identical connectivity) but we're assuming
                    ## that there are not more than one different connectivities
                    ## that have multiple structures
                    #################################################################

                    chosen_connectivity = None
                    for connectivity in connected_structures:
                        if(len(connected_structures[connectivity])>1):
                            #At time of writing, only happens once per compound
                            chosen_connectivity = connectivity

                    #################################################################
                    ## Here, we have a connectivity that is the same for multiple
                    ## structures, so we take those structures and try to pick one
                    ## to use as the primary structure. First we try to see if there 
                    ## is one (and one only) that has stereochemistry. If not,
                    ## again, we prioritize by database
                    #################################################################

                    if(chosen_connectivity is not None):
                        
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

                    #################################################################
                    ## Finally if we get to the point where there are multiple
                    ## connectivities that only have one structure, we prioritize
                    ## by database. There's some redundancy here that can be
                    ## eliminated
                    #################################################################

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

            #################################################################
            ## We have a chosen structure and we collect the aliases for that
            ## structure, so that we can then collect all the structure of
            ## different types for the same aliases
            #################################################################

            chosen_aliases=dict()
            for source in struct_conflicts[chosen_structure]:
                for alias in struct_conflicts[chosen_structure][source]:
                    chosen_aliases[alias]=1

            #################################################################
            ## Having collected aliases for the same structure we iterate
            ## through the types, and the aliases, to find the right structure
            ## to print with the alias
            #################################################################

            for structure_type in "SMILE","InChIKey","InChI":
                if(structure_type not in Structs):
                    continue

                structure_to_use = None
                for alias in chosen_aliases:
                    for structure in Structs[structure_type][struct_stage]:
                        if(alias in Structs[structure_type][struct_stage][structure]):
                            structure_to_use=structure

                #################################################################
                ## Now we collect *all* aliases that are associated with different
                ## structures for the same compound, to associate with the
                ## primary structure
                #################################################################

                aliases=dict()
                for structure in Structs[structure_type][struct_stage]:
                    for alias in Structs[structure_type][struct_stage][structure]:
                        aliases[alias]=1

                #################################################################
                ## Finally, write to file
                #################################################################

                unique_structs_file.write("\t".join((msid,\
                                                     structure_type,\
                                                     ";".join(sorted(aliases)),\
                                                     formula_charge_dict['formula'],\
                                                     formula_charge_dict['charge'],\
                                                     structure_to_use))+"\n")
                                        
    #################################################################
    ## Here we report the structural conflicts
    #################################################################

    if(struct_conflict==1):
        for structure in Structs[struct_type][struct_stage]:
            for external_id in Structs[struct_type][struct_stage][structure]:
                structure_conflicts_file.write("\t".join((msid,struct_type,struct_stage,structure,external_id,
                                                          Structs[struct_type][struct_stage][structure][external_id]))+"\n")
                                        
    #################################################################
    ## Here we report the formula conflicts
    #################################################################

    if(formula_conflict==1):
        for formula in Formulas[struct_type][struct_stage]:
            for external_id in Formulas[struct_type][struct_stage][formula]:
                formula_dict=json.loads(formula)
                formula_conflicts_file.write("\t".join((msid,struct_type,struct_stage,formula_dict['formula'],formula_dict['charge'],external_id,
                                                        Formulas[struct_type][struct_stage][formula][external_id]))+"\n")
