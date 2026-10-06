#!/usr/bin/env python

if __name__ == "__main__":
    # Argument guard -- see "The argument guard" in Scripts/README.md.
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
		cpd = compounds_dict[rgt['compound']]
		# Formula AND charge. This script refreshed only the formula, so a
		# compound whose protonation state changed left its old charge behind
		# in every reaction's stoichiometry; balanceReaction() reads that
		# embedded pair, saw a phantom hydrogen imbalance with no charge
		# imbalance, and Adjust_Reaction_Protons.py "fixed" it by adding a
		# proton -- manufacturing a real charge imbalance on 5,935 reactions
		# that were balanced against the records (OK -> CI:1), and rewriting
		# proton coefficients on 17,457. It never showed before because no
		# earlier refresh recharged 7,500 compounds at once.
		changed = False
		if(rgt['formula'] != cpd['formula']):
			print("Updating formula in stoichiometry for",rxn,"from",rgt['formula'],"to",cpd['formula'])
			rgt['formula'] = cpd['formula']
			changed = True
		if(rgt.get('charge') != cpd['charge']):
			print("Updating charge in stoichiometry for",rxn,"from",rgt.get('charge'),"to",cpd['charge'])
			rgt['charge'] = cpd['charge']
			changed = True
		if(changed and rxn not in updated_reactions_list):
			updated_reactions_list.append(rxn)

if(len(updated_reactions_list)>0):
    print("Saving rebuilt stoichiometries for "+str(len(updated_reactions_list))+" reactions")
    reactions_helper.saveReactions(reactions_dict)