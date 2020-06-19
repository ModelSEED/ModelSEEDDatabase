#!/usr/bin/env bash
date
DIR=`dirname $0`
echo "Rebuild"
${DIR}/Rebuild_Reactions.py
echo "Merge reactions"
${DIR}/Merge_Reactions.py
echo "Rebalance"
${DIR}/Rebalance_Reactions.py
echo "Adjust protons"
${DIR}/Adjust_Reaction_Protons.py
echo "Adjust water"
${DIR}/Adjust_Reaction_Water.py
echo "Merge obsolete aliases"
${DIR}/Merge_Obsolete_Aliases.py
echo "Update aliases"
${DIR}/Update_Compound_Aliases_in_DB.py
${DIR}/Update_Reaction_Aliases_in_DB.py
echo "Refresh complete: use git diff to observe changes"
date