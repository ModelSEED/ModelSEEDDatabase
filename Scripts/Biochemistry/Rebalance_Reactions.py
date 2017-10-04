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

    Rxn_Cpds_Array=ReactionsHelper.parseStoich(Reactions_Dict[rxn]["stoichiometry"])
    Status = ReactionsHelper.balanceReaction(Rxn_Cpds_Array)

    #Need to handle reactions with polymers
    if("ERROR" in Status):
        continue

    old_status=Reactions_Dict[rxn]["status"]
    if(Status != old_status and "CK" not in old_status):
        print "Changing Status for "+rxn+" from "+old_status+" to "+Status
        Reactions_Dict[rxn]["status"]=Status
        Update_Reactions=1

if(Update_Reactions==1):
    print "Saving reactions";
    ReactionsHelper.saveReactions(Reactions_Dict)
