#!/usr/bin/env python

if __name__ == "__main__":
    # Validate before importing or touching the database. 'verbose' was matched
    # by substring against the whole of argv, so any unrecognised argument was
    # silently ignored and the script ran with its defaults.
    import argparse as _argparse
    _p = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter)
    # The body reads three flags out of argv -- "save" (write the new status
    # back; a dry run otherwise), "print" (per-reaction charge detail) and
    # "verbose" -- but the guard admitted only "verbose", so the documented
    # `Rebalance_Reactions.py save` in Refresh_Reactions.sh was rejected with
    # "invalid choice: 'save'" and the rebalance step of the refresh sequence
    # silently never ran. Admit what the body actually consumes.
    _p.add_argument("flags", nargs="*", choices=("verbose", "save", "print"),
                    help="save: write the recomputed status back (dry run "
                         "otherwise); print: per-reaction charge detail; "
                         "verbose: per-reaction detail")
    _p.parse_args()


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

    #Need to handle reactions with polymers
    if(new_status=="Duplicate reagents"):
        new_status = "NB"
        continue

    # CK reactions are rewritten too, not merely warned about. A disagreement
    # on a curator-checked reaction used to print a warning and leave the
    # stale verdict in place, so the field kept asserting OK for reactions
    # that no longer balanced. The CK marker is carried across by preserveCK
    # so the curator record survives the rewrite.
    new_status = reactions_helper.preserveCK(new_status, old_status)

    if(new_status != old_status):
        if("CK" in old_status):
            print("Updating previously checked (CK) reaction "+rxn+": "+old_status+" -> "+new_status)
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
