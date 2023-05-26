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
from BiochemPy import Reactions, Compounds

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

Update_Reactions=0
status_lines = list()
for rxn in sorted(reactions_dict.keys()):
    if(reactions_dict[rxn]["status"] == "EMPTY"):
        continue

    rxn_cpds_array=reactions_dict[rxn]["stoichiometry"]
    
    # Check that all reagents have structures
    all_structures=True
    for rgt in rxn_cpds_array:
        if(compounds_dict[rgt['compound']]['smiles'] == '' and \
            compounds_dict[rgt['compound']]['inchikey'] == ''):
            all_structures=False

    new_status = reactions_helper.balanceReaction(rxn_cpds_array, all_structures)
    old_status=reactions_dict[rxn]["status"]

    if("CK" in old_status and new_status not in old_status):
        print("Warning: previously checked (CK) reaction may have different status: ",rxn,new_status,old_status)
    
    #Need to handle reactions with polymers
    if(new_status=="Duplicate reagents"):
        new_status = "NB"
        continue

    if(new_status != old_status and "CK" not in old_status):
        print("Changing Status for "+rxn+" from "+old_status+" to "+new_status)
        status_lines.append(rxn+"\t"+old_status+"\t"+new_status+"\n")
        reactions_dict[rxn]["status"]=new_status
        Update_Reactions+=1
        if(print_charges is True):
            for entry in rxn_cpds_array:
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
        reactions_helper.saveReactions(reactions_dict)
