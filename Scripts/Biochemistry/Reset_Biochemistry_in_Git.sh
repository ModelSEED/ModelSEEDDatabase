#!/usr/bin/env bash
#Compounds
git checkout ../../Biochemistry/compound_*.tsv
git checkout ../../Biochemistry/compound_*.json
git checkout ../../Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt
git checkout ../../Biochemistry/Aliases/Unique_ModelSEED_Compound_Names.txt

#Reactions
git checkout ../../Biochemistry/reaction_*.tsv
git checkout ../../Biochemistry/reaction_*.json
git checkout ../../Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt
git checkout ../../Biochemistry/Aliases/Unique_ModelSEED_Reaction_Names.txt
git checkout ../../Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt

#Structures
git checkout ../../Biochemistry/Structures/All_ModelSEED_Structures.txt
git checkout ../../Biochemistry/Structures/Unique_ModelSEED_Structures.txt
