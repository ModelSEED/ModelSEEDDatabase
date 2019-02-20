#!/usr/bin/env python
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
