#!/usr/bin/env bash
./Rebuild_Reactions.py
./Merge_Reactions.py
./Rebalance_Reactions.py
./Adjust_Reaction_Protons.py
./Adjust_Reaction_Water.py
./Merge_Obsolete_Aliases.py
./Update_Compound_Aliases_in_DB.py
./Update_Reaction_Aliases_in_DB.py
echo "Refresh complete: use git diff to observe changes"
