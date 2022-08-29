#!/usr/bin/env python
import os
import sys
import subprocess
import time
import copy
import re
import json
from collections import OrderedDict

sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions, Compounds
compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

updated_reactions=list()
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

            rxn_cpds_array = reactions_dict[rxn]["stoichiometry"]
            stoichiometry=reactions_helper.buildStoich(rxn_cpds_array)
            
            reactions_helper.replaceCompound(rxn_cpds_array,cpd,lnkd_cpd)
            new_stoichiometry=reactions_helper.buildStoich(rxn_cpds_array)

            if(stoichiometry != new_stoichiometry):
                print("Replacting obsolete "+cpd+" with "+lnkd_cpd+" in "+rxn)
                updated_reactions.append(rxn)

#    if(Update_Reactions>0):
#        break

if(len(updated_reactions)>0):
    print("Saving replacement of "+str(len(updated_reactions))+" obsolete compounds in reactions")
    reactions_helper.saveReactions(reactions_dict)

# Should probably test for reaction balance, but the merged compounds should have the same structure/formula
# So running Structures/Update_Compound_Structures_Formulas_Charge.py and
# ./Rebalance_Reactions.py would fix any changes, but *should not be needed*

# Having updated a key aspect of the reaction stoichiometry, we must check to see if there's any newly obsolete reactions
# Should run ./Merge_Reactions.py and ./Remove_Newly_Obsolescent_Reactions.py

