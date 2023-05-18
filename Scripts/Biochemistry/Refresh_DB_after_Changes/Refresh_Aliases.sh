#!/usr/bin/env bash
date
DIR=`dirname $0`
echo "Merge obsolete aliases"
${DIR}/Merge_Obsolete_Aliases.py
echo "Update aliases"
${DIR}/Update_Compound_Aliases_in_DB.py
${DIR}/Update_Reaction_Aliases_in_DB.py
echo "Refresh Aliases complete: use git diff to observe changes"
date
