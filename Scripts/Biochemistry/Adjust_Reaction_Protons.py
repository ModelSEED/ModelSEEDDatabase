#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

Update_Reactions=0
for rxn in sorted(Reactions_Dict.keys()):

    if(Reactions_Dict[rxn]["status"] == "EMPTY"):
        continue

    #Find statuses that only have proton imbalance
    if("MI" not in Reactions_Dict[rxn]["status"]):
        continue

    Status_Blocks = Reactions_Dict[rxn]["status"].split("|")
    for block in Status_Blocks:
        if("MI" not in block):
            continue

        block = block.replace("MI:","")
        elements = block.split("/")

        #only making adjustments if mass imbalance is a single element,
        #and the element is hydrogen
        if(len(elements)>1 or not elements[0].startswith("H:")):
            continue

        old_stoichiometry=Reactions_Dict[rxn]["stoichiometry"]
        Rxn_Cpds_Array=ReactionsHelper.parseStoich(old_stoichiometry)
        ReactionsHelper.adjustCompound(Rxn_Cpds_Array,"cpd00067",float(number))
        new_status = ReactionsHelper.balanceReaction(Rxn_Cpds_Array)
        new_stoichiometry = ReactionsHelper.buildStoich(Rxn_Cpds_Array)

        if(new_stoichiometry != old_stoichiometry):
            print "Rebuilding reaction :",rxn
            ReactionsHelper.rebuildReaction(Reactions_Dict[rxn],new_stoichiometry)
            Reactions_Dict[rxn]["status"]=new_status
            if("HB" not in Reactions_Dict[rxn]["notes"]):
                if(Reactions_Dict[rxn]["notes"]=="" or Reactions_Dict[rxn]["notes"]=="null"):
                    Reactions_Dict[rxn]["notes"]="HB"
                else:
                    Reactions_Dict[rxn]["notes"]+="|HB"
            Update_Reactions+=1

if(Update_Reactions>0):
    print "Saving adjusted protons for "+str(Update_Reactions)+" reactions";
    ReactionsHelper.saveReactions(Reactions_Dict)
