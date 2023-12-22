#!/usr/bin/env python
from os import stat
import sys
sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions, Compounds

be_verbose=False
if("verbose" in sys.argv):
    be_verbose=True

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

Update_Reactions=0
status_lines = list()
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
        
        (element,number)=elements[0].split(":")
        if(be_verbose is True):
            print("Adjusting: "+rxn,element,number)

        #Parse old stoichiometry into array
        rgts_array=Reactions_Dict[rxn]["stoichiometry"]
        old_stoichiometry=ReactionsHelper.buildStoich(rgts_array)

        #Adjust for protons
        ReactionsHelper.adjustCompound(rgts_array,"cpd00067",float(number))

        # Check that all reagents have structures
        all_structures=True
        for rgt in rgts_array:
            if(compounds_dict[rgt['compound']]['smiles'] == '' and \
                compounds_dict[rgt['compound']]['inchikey'] == ''):
                all_structures=False

        #Recompute new status and stoichiometry
        new_status = ReactionsHelper.balanceReaction(rgts_array,all_structures)
        new_stoichiometry = ReactionsHelper.buildStoich(rgts_array)

        if(new_status != Reactions_Dict[rxn]['status']):
            status_lines.append(rxn+"\t"+Reactions_Dict[rxn]['status']+"\t"+new_status+"\n")

        if(new_stoichiometry != old_stoichiometry):
            if(be_verbose is True):
                print("Rebuilding reaction :",rxn)
            ReactionsHelper.rebuildReaction(Reactions_Dict[rxn],rgts_array)
            if(be_verbose is True):
                print("Saving new status for reaction :",rxn,new_status)
            Reactions_Dict[rxn]["status"]=new_status
            if("HB" not in Reactions_Dict[rxn]["notes"]):
                Reactions_Dict[rxn]["notes"].append("HB")
            Update_Reactions+=1

if(len(status_lines)>0):
    print("Updating status for "+str(len(status_lines))+" reactions")
    status_file = open("Status_Changes_After_Proton_Adjustment.txt",'w')
    for line in status_lines:
        status_file.write(line)
    status_file.close()

if(Update_Reactions>0):
    print("Saving adjusted protons for "+str(Update_Reactions)+" reactions")
    ReactionsHelper.saveReactions(Reactions_Dict)
