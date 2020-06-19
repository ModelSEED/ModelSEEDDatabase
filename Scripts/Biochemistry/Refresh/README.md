# Order with which these scripts should be run. If nothing has been changed in the
# database, then these scripts will then not change anything, which you can check
# using git diff

# 1) This one attempts to balance "out" compounds, and is key for
# identifying reactions wherein the updated stoichiometry should be
# adjusted for the purposes of FBA.
./Rebuild_Reactions.py

# 2) This one matches reaction stoichiometry and indicates which reactions
# are now obsolete
./Merge_Reactions.py

# 3) This one calculates the new status of the reaction based on how the mass
# and charge balances out in the reactants and products. Its different from
# Rebuild_Reactions.py because it doesn't attempt to change the stoichiometry
# but it does ignore duplicate compounds
./Rebalance_Reactions.py

# 4) These two scripts take the result of the prior script, and uses it to find
# reactions which can be balanced with protons and water. They print out results
# in an output file
./Adjust_Reaction_Protons.py
./Adjust_Reaction_Water.py

# 5) This script takes the output of Merge_Reactions.py and makes sure that all the
# aliases are consistently applied to linked reactions and compounds
./Merge_Obsolete_Aliases.py

# 6) Based on how aliases are updated, they are usually only updated in the 'Aliases' files,
# so these two scripts makes sure that they are updated in the actual database. They have
# some rules for what gets updated so not to overwhelm the user
./Update_Compound_Aliases_in_DB.py
./Update_Reaction_Aliases_in_DB.py

