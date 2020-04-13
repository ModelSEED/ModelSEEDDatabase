#!/usr/bin/env python
from BiochemPy import Reactions, Compounds

##########################################
# Fixing compounds
##########################################
compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
for cpd in compounds_dict:
    break
#compounds_helper.saveCompounds(compounds_dict)

###############################################
# Fixing reactions
###############################################
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
for rxn in reactions_dict:
    break
#reactions_helper.saveReactions(reactions_dict)
