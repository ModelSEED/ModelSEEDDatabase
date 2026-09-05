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
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions, Compounds

##########################################
# Reprinting compounds, aliases and names
##########################################
compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
compounds_helper.saveCompounds(compounds_dict)
aliases_dict = compounds_helper.loadMSAliases()
compounds_helper.saveAliases(aliases_dict)
names_dict = compounds_helper.loadNames()
compounds_helper.saveNames(names_dict)

###############################################
# Reprinting reactions, aliases, names and ecs
###############################################
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
for rxn in reactions_dict:
    reactions_helper.rebuildReaction(reactions_dict[rxn])
reactions_helper.saveReactions(reactions_dict)
aliases_dict = reactions_helper.loadMSAliases()
reactions_helper.saveAliases(aliases_dict)
names_dict = reactions_helper.loadNames()
reactions_helper.saveNames(names_dict)
ecs_dict = reactions_helper.loadECs()
reactions_helper.saveECs(ecs_dict)
