#!/usr/bin/env python
import sys
sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions, Compounds

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
for cpd in compounds_dict:
    # consistently add or alter attribute
    pass
compounds_helper.saveCompounds(compounds_dict)

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
for rxn in reactions_dict:
    # consistently add or alter attribute
    pass
reactions_helper.saveReactions(reactions_dict)