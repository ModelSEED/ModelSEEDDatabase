#!/usr/bin/env python
import sys
sys.path.append('../../../Libs/Python/')
from BiochemPy import Reactions

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

Update_Reactions=0
for rxn in sorted(Reactions_Dict.keys()):
    if(Reactions_Dict[rxn]["status"] == "EMPTY"):
        continue

    Rxn_Cpds_Array = Reactions_Dict[rxn]["stoichiometry"]
    Stoichiometry=ReactionsHelper.buildStoich(Rxn_Cpds_Array)

    New_Rxn_Cpds_Array = ReactionsHelper.removeCpdRedundancy(Rxn_Cpds_Array)
    New_Stoichiometry=ReactionsHelper.buildStoich(New_Rxn_Cpds_Array)

    if(Stoichiometry != New_Stoichiometry):
        print("Rebuilding "+rxn)
        ReactionsHelper.rebuildReaction(Reactions_Dict[rxn], New_Rxn_Cpds_Array)
        Update_Reactions+=1

if(Update_Reactions>0):
    print("Saving rebuilt equations for "+str(Update_Reactions)+" reactions")
    ReactionsHelper.saveReactions(Reactions_Dict)
