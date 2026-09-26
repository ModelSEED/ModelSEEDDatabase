After running Update_Curated_Compound_Attributes.py, you should run these scripts to update the database:

../../Biochemistry/Refresh/Rebalance_Reactions.py (very important)
../../Biochemistry/Refresh/Adjust_Reaction_Protons.py
../../Biochemistry/Refresh/Adjust_Reaction_Water.py
../../Biochemistry/Refresh/Merge_Reactions.py (merges may happen because of water, so need to double-check the results)
../../Biochemistry/Refresh/Update_Compound_Aliases_in_DB.py
