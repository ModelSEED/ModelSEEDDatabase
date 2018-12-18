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
    new_status = ReactionsHelper.balanceReaction(Rxn_Cpds_Array)

    #Need to handle reactions with polymers
    if("ERROR" in new_status):
        continue

    old_status=Reactions_Dict[rxn]["status"]
    if(new_status != old_status and "CK" not in old_status):
        if("OK" not in old_status):
            print rxn+"\t"+old_status+"\t"+new_status
            Reactions_Dict[rxn]["status"]=new_status
            Update_Reactions+=1

if(Update_Reactions>0):
    print "Saving updated statuses for "+str(Update_Reactions)+" reactions";
    ReactionsHelper.saveReactions(Reactions_Dict)
