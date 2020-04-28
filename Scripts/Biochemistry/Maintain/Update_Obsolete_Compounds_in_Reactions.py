#!/usr/bin/env python
import os
import sys
import subprocess
import time
import copy
import re
import json
from collections import OrderedDict
from BiochemPy import Reactions, Compounds

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

Update_Reactions=0
for rxn in reactions_dict:
    if(reactions_dict[rxn]['status']=='EMPTY'):
        continue

    for cpd in reactions_dict[rxn]['compound_ids'].split(';'):
        if(compounds_dict[cpd]['is_obsolete'] == 1):
            if(compounds_dict[cpd]['linked_compound'] == 'null'):
                print("Warning: missing linked compound for obsolete compound: "+cpd)
                continue
            
            lnkd_cpd = sorted(compounds_dict[cpd]['linked_compound'].split(';'))[0]

            # Replace cpd with lnkd_cpd in reaction fields:
            # code, compound_ids, equation, stoichiometry

            old_stoichiometry = reactions_dict[rxn]["stoichiometry"]
            rxn_cpds_array = reactions_helper.parseStoich(old_stoichiometry)
            reactions_helper.replaceCompound(rxn_cpds_array,cpd,lnkd_cpd)
            new_stoichiometry = reactions_helper.buildStoich(rxn_cpds_array)
            reactions_helper.rebuildReaction(reactions_dict[rxn],new_stoichiometry)

            print("Replacting obsolete "+cpd+" with "+lnkd_cpd+" in "+rxn)

            Update_Reactions+=1

#    if(Update_Reactions>0):
#        break

if(Update_Reactions>0):
    print("Saving replacement of "+str(Update_Reactions)+" obsolete compounds in reactions")
    reactions_helper.saveReactions(reactions_dict)

# Should probably test for reaction balance, but the merged compounds should have the same structure/formula
# So running Structures/Update_Compound_Structures_Formulas_Charge.py and
# ./Rebalance_Reactions.py would fix any changes, but *should not be needed*

# Having updated a key aspect of the reaction stoichiometry, we must check to see if there's any newly obsolete reactions
# Should run ./Merge_Reactions.py and ./Remove_Newly_Obsolescent_Reactions.py

