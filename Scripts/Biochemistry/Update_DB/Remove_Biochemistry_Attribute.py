#!/usr/bin/env python
from BiochemPy import Reactions, Compounds
import sys

remove_index=0
remove_string='ontology'

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

for header in range(len(compounds_helper.Headers)):
    if(compounds_helper.Headers[header]==remove_string):
        remove_index=header

del compounds_helper.Headers[remove_index]

for cpd in compounds_dict:
    del compounds_dict[cpd][remove_string]

compounds_helper.saveCompounds(compounds_dict)

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

for header in range(len(reactions_helper.Headers)):
    if(reactions_helper.Headers[header]==remove_string):
        remove_index=header

del reactions_helper.Headers[remove_index]

for rxn in reactions_dict:
    del reactions_dict[rxn][remove_string]

reactions_helper.saveReactions(reactions_dict)
