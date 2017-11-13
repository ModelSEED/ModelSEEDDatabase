#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

Update_Reactions=0
for rxn in sorted(Reactions_Dict.keys()):
    if(Reactions_Dict[rxn]["status"] == "EMPTY"):
        continue

    Rxn_Cpds_Array = ReactionsHelper.parseStoich(Reactions_Dict[rxn]["stoichiometry"])
    Stoichiometry=ReactionsHelper.buildStoich(Rxn_Cpds_Array)
    if(Stoichiometry != Reactions_Dict[rxn]["stoichiometry"]):
        ReactionsHelper.rebuildReaction(Reactions_Dict[rxn],Stoichiometry)
        Update_Reactions=1

print "Saving reactions";
ReactionsHelper.saveReactions(Reactions_Dict)
