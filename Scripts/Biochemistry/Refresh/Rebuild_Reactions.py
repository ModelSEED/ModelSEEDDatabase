#!/usr/bin/env python
from BiochemPy import Reactions

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

Update_Reactions=0
for rxn in sorted(Reactions_Dict.keys()):
    if(Reactions_Dict[rxn]["status"] == "EMPTY"):
        continue

    Rxn_Cpds_Array = ReactionsHelper.parseStoich(Reactions_Dict[rxn]["stoichiometry"])
    New_Rxn_Cpds_Array = ReactionsHelper.removeCpdRedundancy(Rxn_Cpds_Array)
    Stoichiometry=ReactionsHelper.buildStoich(New_Rxn_Cpds_Array)
    if(Stoichiometry != Reactions_Dict[rxn]["stoichiometry"]):
        print("Rebuilding "+rxn)
        ReactionsHelper.rebuildReaction(Reactions_Dict[rxn],Stoichiometry)
        Update_Reactions+=1

if(Update_Reactions>0):
    print("Saving rebuilt equations for "+str(Update_Reactions)+" reactions")
    ReactionsHelper.saveReactions(Reactions_Dict)
