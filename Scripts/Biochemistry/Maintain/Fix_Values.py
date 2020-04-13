#!/usr/bin/env python
from BiochemPy import Reactions, Compounds

##########################################
# Fixing compounds
##########################################
#compounds_helper = Compounds()
#compounds_dict = compounds_helper.loadCompounds()
#for cpd in compounds_dict:
#    if(compounds_dict[cpd]['linked_compound']==""):
#        compounds_dict[cpd]['linked_compound']=None
compounds_helper.saveCompounds(compounds_dict)

###############################################
# Fixing reactions
###############################################
#reactions_helper = Reactions()
#reactions_dict = reactions_helper.loadReactions()

#rxn = 'rxn42111'
#notes = ['HB','GFP','EQP']
#reactions_dict[rxn]['notes'] = notes

#for rxn in reactions_dict:
#    if(reactions_dict[rxn]['deltag'] == 10000000):
#        reactions_dict[rxn]['deltag'] = 'null'
#    if(reactions_dict[rxn]['deltagerr'] == 10000000):
#        reactions_dict[rxn]['deltagerr'] = 'null'

#reactions_helper.saveReactions(reactions_dict)
