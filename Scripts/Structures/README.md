# This script takes all the InChI and SMILE strings from the original sources and extracts the formulas and charges, including protonation states
./Print_Structure_Formula_Charge.py

# This script takes the list of InChI and SMILE strings from the original sources and attempts to consolidate them for the ModelSEED database
# Because it checks for formula conflicts, it needs the output of the previous script
# and it reports structure and formula conflicts in files in the same directory
# Its key output files are ../../Biochemistry/Structures/All_ModelSEED_Structures.txt
# and ../../Biochemistry/Structures/Unique_ModelSEED_Structures.txt
# The latter file is the structural heart, and our goal is to expand on it via curation of the conflicts and integration of more structures
./List_ModelSEED_Structures.py

# This script takes the output of the previous two scripts, and uses them to update the ModelSEED database
./Update_Compound_Structures_Formulas_Charge.py

# These two scripts build tab-separated files of conflicts and redundancies for reporting
# Their output isn't necessarily optimal, or include all possible problems, but is used as a starting point
# for loading into spreadsheets with the goal of curating by a team
# Redundant structures are structures that are found assigned to more than one compound, and could be merged
# Concflicting structures are structures found assigned to a single compound and could be disambiguated
./Report_Redundant_Structures.py
./Report_Conflicting_Structures.py
