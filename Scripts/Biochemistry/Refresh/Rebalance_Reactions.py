#!/usr/bin/env python
import os, sys
temp=list();
header=1;

dry_run = True
if("save" in sys.argv):
    dry_run = False

print_charges = False
if("print" in sys.argv):
    print_charges = True

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

Update_Reactions=0
status_file = open("Status_Changes.txt",'w')
for rxn in sorted(Reactions_Dict.keys()):
    if(Reactions_Dict[rxn]["status"] == "EMPTY"):
        continue

    Rxn_Cpds_Array=ReactionsHelper.parseStoich(Reactions_Dict[rxn]["stoichiometry"])
    new_status = ReactionsHelper.balanceReaction(Rxn_Cpds_Array)
    old_status=Reactions_Dict[rxn]["status"]

    #Need to handle reactions with polymers
    if(new_status=="Duplicate reagents"):
        new_status = "NB"
        continue

    if(new_status != old_status and "CK" not in old_status):
        print("Changing Status for "+rxn+" from "+old_status+" to "+new_status)
        status_file.write(rxn+"\t"+old_status+"\t"+new_status+"\n")
        Reactions_Dict[rxn]["status"]=new_status
        Update_Reactions+=1
        if(print_charges is True):
            for entry in Rxn_Cpds_Array:
                print("\t".join(["\t",entry['compound'],str(entry['coefficient']), \
                                     str(entry['charge']),entry['name']]))

if(Update_Reactions>0):
    print("Updating statuses for "+str(Update_Reactions)+" reactions")
    if(dry_run is False):
        ReactionsHelper.saveReactions(Reactions_Dict)
status_file.close()
