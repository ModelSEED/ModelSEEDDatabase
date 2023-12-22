#!/usr/bin/env python
from BiochemPy import Reactions, Compounds
import sys

##########################################
# Reprinting compounds, aliases and names
##########################################
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
#reactions_helper.Headers.append('ontology')
for rxn in reactions_dict:
#    reactions_dict[rxn]['ontology']='class:null|context:null|step:null'
#    reactions_dict[rxn]['ontology']="null"
    if('None' in reactions_dict[rxn]['notes']):
        reactions_dict[rxn]['notes']="null"
reactions_helper.saveReactions(reactions_dict)

#compounds_helper = Compounds()
#compounds_dict = compounds_helper.loadCompounds()
#compounds_helper.Headers.append('notes')
#for cpd in compounds_dict:
#    compounds_dict[cpd]['notes']="null"
#compounds_helper.saveCompounds(compounds_dict)
