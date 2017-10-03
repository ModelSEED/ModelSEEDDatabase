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

Compound_To_Merge_From="cpd19013" #Ammonium
Compound_To_Merge_To="cpd00013" #Ammonia

Update_Reactions=0
for rxn in sorted(Reactions_Dict.keys()):
    if(Reactions_Dict[rxn]["status"] == "EMPTY"):
        continue

    if(Reactions_Dict[rxn]["status"] == "CPDFORMERROR"):
        continue

    old_stoichiometry=Reactions_Dict[rxn]["stoichiometry"]
    Rxn_Cpds_Array=ReactionsHelper.parseStoich(old_stoichiometry)

    to_replace=0
    for entry in Rxn_Cpds_Array:
        if(Compound_To_Merge_From in entry['compound']):
            to_replace=entry['coefficient']

    if(to_replace==0):
        continue

    #Remove old compound, it deletes itself if its coefficient becomes zero
    ReactionsHelper.adjustCompound(Rxn_Cpds_Array,Compound_To_Merge_From,to_replace)

    #Add new compound with same coefficient
    ReactionsHelper.adjustCompound(Rxn_Cpds_Array,Compound_To_Merge_To,to_replace)

    new_status = ReactionsHelper.balanceReaction(Rxn_Cpds_Array)
    new_stoichiometry = ReactionsHelper.buildStoich(Rxn_Cpds_Array)

    print "Rebuilding reaction :",rxn
    ReactionsHelper.rebuildReaction(Reactions_Dict[rxn],new_stoichiometry)
    Reactions_Dict[rxn]["status"]=new_status
    if("ME" not in Reactions_Dict[rxn]["notes"]):
        if(Reactions_Dict[rxn]["notes"]=="" or Reactions_Dict[rxn]["notes"]=="null"):
            Reactions_Dict[rxn]["notes"]="ME"
        else:
            Reactions_Dict[rxn]["notes"]+="|ME"
    Update_Reactions=1
            
    break

#4) You need to update Aliases
#5) You need to update media
#6) You need to report possible updates in templates (and, following modifications, re-build any public models?)

if(Update_Reactions==1):
    print "Saving reactions";
    ReactionsHelper.saveReactions(Reactions_Dict)
