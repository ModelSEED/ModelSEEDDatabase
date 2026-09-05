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
from BiochemPy import Reactions, Compounds

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

updated_reactions=list()
for rxn in reactions_dict:
    if(reactions_dict[rxn]["status"] == "EMPTY"):
        continue

    rxn_cpds_array = reactions_dict[rxn]["stoichiometry"]
    stoichiometry=reactions_helper.buildStoich(rxn_cpds_array)

    # Before we start removing compounds we want to check that the 
    # list of reagents has structural integrity, so that all
    # atoms could otherwise be mapped despite any redundancy
    all_structures=True
    for rgt in rxn_cpds_array:
        if(compounds_dict[rgt['compound']]['smiles'] == '' and \
            compounds_dict[rgt['compound']]['inchikey'] == ''):
            all_structures=False

    if(all_structures is False):
        new_rxn_cpds_array = reactions_helper.removeCpdRedundancy(rxn_cpds_array)
        new_stoichiometry=reactions_helper.buildStoich(new_rxn_cpds_array)

        if(stoichiometry != new_stoichiometry):
            print("Rebuilding "+rxn)
            reactions_helper.rebuildReaction(reactions_dict[rxn], new_rxn_cpds_array)
            updated_reactions.append(rxn)

if(len(updated_reactions)>0):
    print("Saving rebuilt equations for "+str(len(updated_reactions))+" reactions")
    reactions_helper.saveReactions(reactions_dict)
