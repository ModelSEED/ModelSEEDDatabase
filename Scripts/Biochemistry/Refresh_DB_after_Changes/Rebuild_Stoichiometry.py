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

updated_reactions_list=list()
for rxn in reactions_dict:
	if(reactions_dict[rxn]["status"] == "EMPTY"):
		continue
	
	for rgt in reactions_dict[rxn]['stoichiometry']:
		if(rgt['formula'] != compounds_dict[rgt['compound']]['formula']):
			print("Updating formula in stoichiometry for",rxn,"from",rgt['formula'],"to",compounds_dict[rgt['compound']]['formula'])
			rgt['formula'] = compounds_dict[rgt['compound']]['formula']
			if(rxn not in updated_reactions_list):
				updated_reactions_list.append(rxn)

if(len(updated_reactions_list)>0):
    print("Saving rebuilt stoichiometries for "+str(len(updated_reactions_list))+" reactions")
    reactions_helper.saveReactions(reactions_dict)