#!/usr/bin/env python

if __name__ == "__main__":
    # Validate arguments BEFORE importing anything or touching the database.
    # These scripts mutate the database, and without this an unknown flag or a
    # mistyped mode was silently ignored and the script ran with its defaults:
    # asking Estimate_Reaction_Reversibility.py for --help rewrote 122 files.
    # Placed above the imports so --help works even where a dependency is
    # missing from the path.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import sys
sys.path.append('../../../Libs/Python/')
from BiochemPy import Reactions

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

Update_Reactions=0
for rxn in sorted(Reactions_Dict.keys()):

    if(Reactions_Dict[rxn]["status"] == "EMPTY"):
        continue

    old_stoichiometry=Reactions_Dict[rxn]["stoichiometry"]
    is_transport = ReactionsHelper.isTransport(old_stoichiometry)

    if(is_transport != Reactions_Dict[rxn]["is_transport"]):
        print("Updating: ",rxn,is_transport)
        Reactions_Dict[rxn]["is_transport"] = is_transport
        Update_Reactions+=1

if(Update_Reactions>0):
    print("Saving adjusted is_transport flag for "+str(Update_Reactions)+" reactions")
    ReactionsHelper.saveReactions(Reactions_Dict)
