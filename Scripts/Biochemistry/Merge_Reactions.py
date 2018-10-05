#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()
Reactions_Codes = ReactionsHelper.generateCodes(Reactions_Dict)

Update_Reactions=0
for code in sorted(Reactions_Codes.keys()):
    
    primary_rxn = sorted(Reactions_Codes[code].keys())[0]

    if(len(Reactions_Codes[code].keys())==1):
        if(Reactions_Dict[primary_rxn]["is_obsolete"]==1):
            Update_Reactions+=1
            Reactions_Dict[primary_rxn]["is_obsolete"]=0
            Reactions_Dict[primary_rxn]["linked_reaction"]="null"
    else:
        for rxn in Reactions_Codes[code].keys():
            rxn_list = ";".join(sorted(x for x in Reactions_Codes[code].keys() if x != rxn))
            if(rxn == primary_rxn and Reactions_Dict[rxn]["is_obsolete"]==1):
                print "Primary: ",rxn,code,rxn_list,Reactions_Dict[rxn]["linked_reaction"],Reactions_Dict[rxn]["is_obsolete"]
                Update_Reactions+=1
                Reactions_Dict[rxn]["linked_reaction"] = rxn_list
                Reactions_Dict[rxn]["is_obsolete"]=0
            elif(rxn != primary_rxn and Reactions_Dict[rxn]["is_obsolete"]==0):
                print "Secondary: ",rxn,code,rxn_list,Reactions_Dict[rxn]["linked_reaction"],Reactions_Dict[rxn]["is_obsolete"]
                Update_Reactions+=1
                Reactions_Dict[rxn]["linked_reaction"] = rxn_list
                Reactions_Dict[rxn]["is_obsolete"]=1

if(Update_Reactions>0):
    ReactionsHelper.saveReactions(Reactions_Dict)
    pass

