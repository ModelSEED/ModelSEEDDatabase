#!/usr/bin/env python
import sys
sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions, Compounds

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

Update_Reactions=0
status_lines = list()
for rxn in sorted(Reactions_Dict.keys()):
    #Find statuses that only have water imbalance
    if("MI:H:2/O:1" != Reactions_Dict[rxn]["status"] and
       "MI:H:-2/O:-1" != Reactions_Dict[rxn]["status"]):
        continue

    #Parse old stoichiometry into array
    rgts_array=Reactions_Dict[rxn]["stoichiometry"]
    old_stoichiometry=ReactionsHelper.buildStoich(rgts_array)

    #Don't adjust reactions that only have water
    if(len(rgts_array)==1):
        continue

    Water_Adjustment = 1
    if("-1" in Reactions_Dict[rxn]["status"]):
        Water_Adjustment = -1

    #Adjust for water
    ReactionsHelper.adjustCompound(rgts_array,"cpd00001",float(Water_Adjustment))

    # Check that all reagents have structures
    all_structures=True
    for rgt in rgts_array:
        if(compounds_dict[rgt['compound']]['smiles'] == '' and \
            compounds_dict[rgt['compound']]['inchikey'] == ''):
            all_structures=False

    #Recompute new status and stoichiometry
    new_status = ReactionsHelper.balanceReaction(rgts_array, all_structures)
    new_stoichiometry = ReactionsHelper.buildStoich(rgts_array)

    if(new_status != Reactions_Dict[rxn]['status']):
        status_lines.append(rxn+"\t"+Reactions_Dict[rxn]['status']+"\t"+new_status+"\n")

    if(new_stoichiometry != old_stoichiometry):
        print("Rebuilding reaction :",rxn)
        ReactionsHelper.rebuildReaction(Reactions_Dict[rxn],rgts_array)
        Reactions_Dict[rxn]["status"]=new_status
        if("WB" not in Reactions_Dict[rxn]["notes"]):
            if(Reactions_Dict[rxn]["notes"]=="" or Reactions_Dict[rxn]["notes"]=="null"):
                Reactions_Dict[rxn]["notes"]="WB"
            else:
                Reactions_Dict[rxn]["notes"]+="|WB"
        Update_Reactions+=1

if(len(status_lines)>0):
    print("Updating status for "+str(len(status_lines))+" reactions")
    status_file = open("Status_Changes_After_Water_Adjustment.txt",'w')
    for line in status_lines:
        status_file.write(line)
    status_file.close()

if(Update_Reactions>0):
    print("Saving adjusted water for "+str(Update_Reactions)+" reactions")
    ReactionsHelper.saveReactions(Reactions_Dict)
