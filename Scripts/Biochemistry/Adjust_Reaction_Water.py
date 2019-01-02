#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

Update_Reactions=0
status_file = open("Status_Changes_After_Water_Adjustment.txt",'w')
for rxn in sorted(Reactions_Dict.keys()):
    #Find statuses that only have water imbalance
    if("MI:H:2/O:1" != Reactions_Dict[rxn]["status"] and
       "MI:H:-2/O:-1" != Reactions_Dict[rxn]["status"]):
        continue

    #Parse old stoichiometry into array
    old_stoichiometry=Reactions_Dict[rxn]["stoichiometry"]
    Rxn_Cpds_Array=ReactionsHelper.parseStoich(old_stoichiometry)

    #Don't adjust reactions that only have water
    if(len(Rxn_Cpds_Array)==1):
        continue

    Water_Adjustment = 1
    if("-1" in Reactions_Dict[rxn]["status"]):
        Water_Adjustment = -1

    #Adjust for water
    ReactionsHelper.adjustCompound(Rxn_Cpds_Array,"cpd00001",float(Water_Adjustment))

    #Recompute new status and stoichiometry
    new_status = ReactionsHelper.balanceReaction(Rxn_Cpds_Array)
    new_stoichiometry = ReactionsHelper.buildStoich(Rxn_Cpds_Array)

    if(new_status != Reactions_Dict[rxn]['status']):
        status_file.write(rxn+"\t"+Reactions_Dict[rxn]['status']+"\t"+new_status+"\n")

    if(new_stoichiometry != old_stoichiometry):
        print("Rebuilding reaction :",rxn)
        ReactionsHelper.rebuildReaction(Reactions_Dict[rxn],new_stoichiometry)
        Reactions_Dict[rxn]["status"]=new_status
        if("WB" not in Reactions_Dict[rxn]["notes"]):
            if(Reactions_Dict[rxn]["notes"]=="" or Reactions_Dict[rxn]["notes"]=="null"):
                Reactions_Dict[rxn]["notes"]="WB"
            else:
                Reactions_Dict[rxn]["notes"]+="|WB"
        Update_Reactions+=1

if(Update_Reactions>0):
    print("Saving adjusted water for "+str(Update_Reactions)+" reactions")
    ReactionsHelper.saveReactions(Reactions_Dict)
