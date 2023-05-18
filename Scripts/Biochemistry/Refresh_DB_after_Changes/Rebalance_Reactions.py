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

sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

Update_Reactions=0
status_lines = list()
for rxn in sorted(Reactions_Dict.keys()):
    if(Reactions_Dict[rxn]["status"] == "EMPTY"):
        continue

    Rxn_Cpds_Array=Reactions_Dict[rxn]["stoichiometry"]
    new_status = ReactionsHelper.balanceReaction(Rxn_Cpds_Array)
    old_status=Reactions_Dict[rxn]["status"]

    #Need to handle reactions with polymers
    if(new_status=="Duplicate reagents"):
        new_status = "NB"
        continue

    if(new_status != old_status and "CK" not in old_status):
        print("Changing Status for "+rxn+" from "+old_status+" to "+new_status)
        status_lines.append(rxn+"\t"+old_status+"\t"+new_status+"\n")
        Reactions_Dict[rxn]["status"]=new_status
        Update_Reactions+=1
        if(print_charges is True):
            for entry in Rxn_Cpds_Array:
                print("\t".join(["\t",entry['compound'],str(entry['coefficient']), \
                                     str(entry['charge']),entry['name']]))

if(len(status_lines)>0):
    print("Updating status for "+str(len(status_lines))+" reactions")
    status_file = open("Status_Changes.txt",'w')
    for line in status_lines:
        status_file.write(line)
    status_file.close()

if(Update_Reactions>0):
    print("Updating statuses for "+str(Update_Reactions)+" reactions")
    if(dry_run is False):
        ReactionsHelper.saveReactions(Reactions_Dict)
