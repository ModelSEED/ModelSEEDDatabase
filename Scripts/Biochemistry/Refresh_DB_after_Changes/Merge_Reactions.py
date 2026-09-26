#!/usr/bin/env python

if __name__ == "__main__":
    # Argument guard -- see "The argument guard" in Scripts/README.md.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import sys
sys.path.append('../../../Libs/Python/')
from BiochemPy import Reactions

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
            if(rxn == primary_rxn and ( Reactions_Dict[rxn]["is_obsolete"]==1 or rxn_list != Reactions_Dict[rxn]["linked_reaction"] )):
                print("Updating primary reaction "+rxn+" and removing any indication that its obsolete")
                Update_Reactions+=1
                Reactions_Dict[rxn]["linked_reaction"] = rxn_list
                Reactions_Dict[rxn]["is_obsolete"]=0
            elif(rxn != primary_rxn and ( Reactions_Dict[rxn]["is_obsolete"]==0 or rxn_list != Reactions_Dict[rxn]["linked_reaction"] )):
                print("Updating reaction "+rxn+" to indicate that its obsolete")
                Update_Reactions+=1
                Reactions_Dict[rxn]["linked_reaction"] = rxn_list
                Reactions_Dict[rxn]["is_obsolete"]=1

if(Update_Reactions>0):
    print("Saving obsolescence updating in "+str(Update_Reactions)+" reactions")
    print("Run ./Merge_Obsolete_Aliases.py and ./Update_Reaction_Aliases.py")
    ReactionsHelper.saveReactions(Reactions_Dict)



