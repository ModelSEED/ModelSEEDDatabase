#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()

Compound_To_Merge_From="cpd00013"
Compound_To_Merge_To="cpd19013"

Cpds_Rxns_Dict=dict()
Rxns_Cpds_Dict=dict()
for rxn in Reactions_Dict.keys():
    if(Reactions_Dict[rxn]["status"] == "EMPTY"):
        continue

    for rgt in Reactions_Dict[rxn]["stoichiometry"].split(";"):
        (coeff,cpd,cpt,index,name) = rgt.split(":",4)

        if(cpd not in Cpds_Rxns_Dict):
            Cpds_Rxns_Dict[cpd]=dict()
        Cpds_Rxns_Dict[cpd][rxn]=1

        if(rxn not in Rxns_Cpds_Dict):
            Rxns_Cpds_Dict[rxn]=dict()
        Rxns_Cpds_Dict[rxn][cpd]=1

#Merging two compounds means:
#1) You take all the reactions for the second compound, and replace the compound id in the second reaction with the first compound
        #Change stoichiometry only, first
#2) You check to see if all the reactions are balanced following the change
#3) You need to check to see if new reactions are now merged/linked to other reactions in database
        #If truly still new reaction, change definition, code, compound_ids, equation
        #If merged, change to obsolete and store linked reaction(s)
#4) You need to update Aliases
#5) You need to update media
#6) You need to report possible updates in templates (and, following modifications, re-build any public models?)
 
for rxn in Cpds_Rxns_Dict[Compound_To_Merge_From].keys():
    old_stoichiometry = Reactions_Dict[rxn]["stoichiometry"]
    new_stoichiometry_array = list()
    for rgt in old_stoichiometry.split(";"):
        (coeff,cpd,cpt,index,name) = rgt.split(":",4)

        #Replace cpd
        if(cpd == Compound_To_Merge_From):
            cpd = Compound_To_Merge_To
        new_stoichiometry_array.append(":".join([coeff,cpd,cpt,index,name]))
    new_stoichiometry = ";".join(new_stoichiometry_array)

    if(new_stoichiometry == old_stoichiometry):
        print rxn, old_stoichiometry, new_stoichiometry
    break



sys.exit()
Update_Reactions=0
for rxn in sorted(Reactions_Dict.keys()):
    if(Reactions_Dict[rxn]["status"] == "EMPTY"):
        continue

    Rxn_Cpds_Array=list()
    for rgt in Reactions_Dict[rxn]["stoichiometry"].split(";"):
        (coeff,cpd,cpt,index,name) = rgt.split(":",4)
        rgt_id = cpd+"_"+cpt+index
        Rxn_Cpds_Array.append({"reagent":rgt_id,"coefficient":coeff,
                               "formula":Compounds_Dict[cpd]["formula"],
                               "charge":Compounds_Dict[cpd]["charge"]})
    Status = ReactionsHelper.balanceReaction(Rxn_Cpds_Array)

    if("ERROR" in Status):
#        print rxn,Status
        continue

    #Remove old HB message
    old_status = ""
    for item in Reactions_Dict[rxn]["status"].split("|"):
        if(item != "HB"):
            old_status += item+"|"
    old_status = old_status[0:-1]

    if(Status != old_status and ("OK" not in old_status and "OK" not in Status)):
        print "Changing Status for "+rxn+" from "+Reactions_Dict[rxn]["status"]+" to "+Status
        Reactions_Dict[rxn]["status"]=Status
        Update_Reactions=1

#if(Update_Reactions==1):
#    print "Saving reactions";
#    ReactionsHelper.saveReactions(Reactions_Dict)
