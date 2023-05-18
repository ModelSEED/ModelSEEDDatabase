#!/usr/bin/env python
import sys
sys.path.append('../../../Libs/Python/')
from BiochemPy import Reactions

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

updated_reactions=list()
for rxn in sorted(reactions_dict.keys()):
    if(reactions_dict[rxn]["status"] == "EMPTY"):
        continue

    rxn_cpds_array = reactions_dict[rxn]["stoichiometry"]
    stoichiometry=reactions_helper.buildStoich(rxn_cpds_array)

    new_rxn_cpds_array = reactions_helper.removeCpdRedundancy(rxn_cpds_array)
    new_stoichiometry=reactions_helper.buildStoich(new_rxn_cpds_array)

    if(stoichiometry != new_stoichiometry):
        print("Rebuilding "+rxn)
        reactions_helper.rebuildReaction(reactions_dict[rxn], new_rxn_cpds_array)
        updated_reactions.append(rxn)

if(len(updated_reactions)>0):
    print("Saving rebuilt equations for "+str(len(updated_reactions))+" reactions")
    reactions_helper.saveReactions(reactions_dict)
