First, you run Add_New_Curated_Compounds.py, it takes two positional
arguments:
1) The tab-separated flat-file of compounds with their names,
structures, and other attributes
2) The GitHub username for the curator adding the new compounds

If you use the '-r' flag, it generates a simple report of what matches
and why, you can use it to modify your file so that you get better
matches. In order to run the script to add new reactions, you have to
_save_ the new compounds in the database using the '-s' flag.

If you do save, but find you have to revert, you'll have to use git to
checkout everything in the Biochemistry folder, which is why it's
recommended to generate the report first

An example:
./Add_New_Curated_Compounds.py ../../../Biochemistry/Aliases/Provenance/ChEBI_ID_Name_InChIKey.tsv samseaver -r -s

Then you run Add_New_Curated_Reactions.py, it takes three positional
arguments:

1) The tab-separated flat-file of reactions with their equations
2) the name of the user or external database within the biochemistry
to use to match the compounds, this can be different from the source
of the reactions themselves (for example Rhea uses ChEBI compound
identifiers for their reactions)
3) The GitHub username for the curator adding the new reactions. This
can be the same entry as the second argument

You use the '-r' and the '-s' flags in the same manner to iteratively
improve your matches and then save to the biochemistry

An example:
./Add_New_Curated_Reactions.py ../../../Biochemistry/Aliases/Provenance/Rhea_reactions.tsv ChEBI samseaver -r -s

After running Add_New_Curated_Compounds.py, you should run these scripts to update the database:

../../Biochemistry/Update/Merge_Formulas.py
../../Biochemistry/Refresh/Update_Compound_Aliases_in_DB.py
../../Structures/List_ModelSEED_Structures.py
../../Structures/Update_Compound_Structures_Formulas_Charge.py

After running Add_New_Curated_Reactions.py, you should run these scripts to update the database:

../../Biochemistry/Refresh/Rebalance_Reactions.py (very important)
../../Biochemistry/Refresh/Adjust_Reaction_Protons.py
../../Biochemistry/Refresh/Adjust_Reaction_Water.py
../../Biochemistry/Refresh/Merge_Reactions.py (merges may happen because of water, so need to double-check the results)
../../Biochemistry/Refresh/Update_Reaction_Aliases.py
../../Biochemistry/Refresh/Update_Reaction_Aliases_in_DB.py
