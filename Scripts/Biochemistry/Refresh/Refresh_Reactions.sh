#!/usr/bin/env bash
date
DIR=`dirname $0`
echo "Rebuild"
${DIR}/Rebuild_Reactions.py
echo "Merge reactions"
${DIR}/Merge_Reactions.py
echo "Rebalance"
${DIR}/Rebalance_Reactions.py save
echo "Adjust protons"
${DIR}/Adjust_Reaction_Protons.py
echo "Adjust water"
${DIR}/Adjust_Reaction_Water.py
echo "Refresh Reactions complete: use git diff to observe changes"
date
