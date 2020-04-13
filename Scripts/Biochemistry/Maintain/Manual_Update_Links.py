#!/usr/bin/env python
import os
import sys
import subprocess
import time
import copy
import re
import json
from collections import OrderedDict
from BiochemPy import Reactions, Compounds

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
cpds_aliases_dict = compounds_helper.loadMSAliases()
cpds_names_dict = compounds_helper.loadNames()
structures_dict = compounds_helper.loadStructures(["InChI","SMILE"],["ModelSEED"])

# cpd00000 shouldn't be in aliases_dict
# del(cpds_aliases_dict['cpd00000'])
# compounds_helper.saveAliases(cpds_aliases_dict)

# Load Reactions
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
rxns_aliases_dict = reactions_helper.loadMSAliases()
rxns_names_dict = reactions_helper.loadNames()
rxns_ecs_dict = reactions_helper.loadECs()

"""
rxn41515
         rxn18357 MetaCyc 3.6.3.17-RXN
         rxn18357 metanetx.reaction MNXR115558
"""
rxn='rxn18357'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('3.6.3.17-RXN')
if('metanetx.reaction' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['metanetx.reaction']=list()
rxns_aliases_dict[rxn]['metanetx.reaction'].append('MNXR115558')

del(rxns_aliases_dict['rxn41515'])
"""
rxn41648
         rxn18664 MetaCyc 3.6.3.31-RXN
         rxn18664 metanetx.reaction MNXR115564
         rxn18664 7.6.2.k
"""
rxn='rxn18664'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('3.6.3.31-RXN')
if('metanetx.reaction' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['metanetx.reaction']=list()
rxns_aliases_dict[rxn]['metanetx.reaction'].append('MNXR115564')

# if(rxn not in rxns_ecs_dict):
#    rxns_ecs_dict[rxn]=list()
# rxns_ecs_dict[rxn].append('7.6.2.k')

del(rxns_aliases_dict['rxn41648'])
del(rxns_ecs_dict['rxn41648'])
"""
rxn41786
         rxn18580 MetaCyc 3.6.3.22-RXN
         rxn18580 metanetx.reaction MNXR115561
         rxn18580 7.4.2.2
"""
rxn='rxn18580'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('3.6.3.22-RXN')
if('metanetx.reaction' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['metanetx.reaction']=list()
rxns_aliases_dict[rxn]['metanetx.reaction'].append('MNXR115561')

# if(rxn not in rxns_ecs_dict):
#    rxns_ecs_dict[rxn]=list()
# rxns_ecs_dict[rxn].append('7.4.2.2')

del(rxns_aliases_dict['rxn41786'])
del(rxns_ecs_dict['rxn41786'])
"""
rxn43490
         rxn18480 MetaCyc 3.6.3.18-RXN
         rxn18480 metanetx.reaction MNXR115559
"""
rxn='rxn18480'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('3.6.3.18-RXN')
if('metanetx.reaction' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['metanetx.reaction']=list()
rxns_aliases_dict[rxn]['metanetx.reaction'].append('MNXR115559')

del(rxns_aliases_dict['rxn43490'])
"""
rxn43764
         rxn27627 MetaCyc 3.6.3.18-RXN
         rxn27627 metanetx.reaction MNXR115652
         rxn27627 7.2.2.p
"""
rxn='rxn27627'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('3.6.3.18-RXN')
if('metanetx.reaction' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['metanetx.reaction']=list()
rxns_aliases_dict[rxn]['metanetx.reaction'].append('MNXR115652')

# if(rxn not in rxns_ecs_dict):
#    rxns_ecs_dict[rxn]=list()
# rxns_ecs_dict[rxn].append('7.2.2.p')

del(rxns_aliases_dict['rxn43764'])
del(rxns_ecs_dict['rxn43764'])
"""
rxn43962
         rxn18697 MetaCyc 3.6.3.44-RXN
         rxn18697 metanetx.reaction MNXR115574
"""
rxn='rxn18697'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('3.6.3.44-RXN')
if('metanetx.reaction' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['metanetx.reaction']=list()
rxns_aliases_dict[rxn]['metanetx.reaction'].append('MNXR115574')

del(rxns_aliases_dict['rxn43962'])
"""
rxn45108
         rxn18576 MetaCyc 3.6.3.21-RXN
         rxn18576 metanetx.reaction MNXR115560
"""
rxn='rxn18576'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('3.6.3.21-RXN')
if('metanetx.reaction' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['metanetx.reaction']=list()
rxns_aliases_dict[rxn]['metanetx.reaction'].append('MNXR115560')

del(rxns_aliases_dict['rxn45108'])
"""
rxn45131
         rxn18690 MetaCyc 3.6.3.42-RXN
         rxn18690 metanetx.reaction MNXR115572
         rxn18690 7.5.2.d
"""
rxn='rxn18690'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('3.6.3.42-RXN')
if('metanetx.reaction' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['metanetx.reaction']=list()
rxns_aliases_dict[rxn]['metanetx.reaction'].append('MNXR115572')

# if(rxn not in rxns_ecs_dict):
#    rxns_ecs_dict[rxn]=list()
# rxns_ecs_dict[rxn].append('7.5.2.d')

del(rxns_aliases_dict['rxn45131'])
del(rxns_ecs_dict['rxn45131'])
"""
rxn45335
         rxn27371 MetaCyc RXN-15598
         rxn27371 metanetx.reaction MNXR140566
         rxn27371 7.2.1.d
"""
rxn='rxn27371'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('RXN-15598')
if('metanetx.reaction' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['metanetx.reaction']=list()
rxns_aliases_dict[rxn]['metanetx.reaction'].append('MNXR140566')

# if(rxn not in rxns_ecs_dict):
#    rxns_ecs_dict[rxn]=list()
# rxns_ecs_dict[rxn].append('7.2.1.d')

del(rxns_aliases_dict['rxn45335'])
del(rxns_ecs_dict['rxn45335'])
"""
rxn47330
         rxn18679 MetaCyc 3.6.3.38-RXN
         rxn18679 metanetx.reaction MNXR115568
         rxn18679 7.6.2.n
"""
rxn='rxn18679'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('3.6.3.38-RXN')
if('metanetx.reaction' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['metanetx.reaction']=list()
rxns_aliases_dict[rxn]['metanetx.reaction'].append('MNXR115568')

# if(rxn not in rxns_ecs_dict):
#    rxns_ecs_dict[rxn]=list()
# rxns_ecs_dict[rxn].append('7.6.2.n')

del(rxns_aliases_dict['rxn47330'])
del(rxns_ecs_dict['rxn47330'])
"""
rxn48286
         rxn18907 MetaCyc 3.6.3.48-RXN
         rxn18907 metanetx.reaction MNXR115576
         rxn18907 7.4.2.g
"""
rxn='rxn18907'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('3.6.3.48-RXN')
if('metanetx.reaction' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['metanetx.reaction']=list()
rxns_aliases_dict[rxn]['metanetx.reaction'].append('MNXR115576')

# if(rxn not in rxns_ecs_dict):
#    rxns_ecs_dict[rxn]=list()
# rxns_ecs_dict[rxn].append('7.4.2.g')

del(rxns_aliases_dict['rxn48286'])
del(rxns_ecs_dict['rxn48286'])
"""
rxn42640
         rxn02157 MetaCyc RXN-115
         rxn02157 rhea 15006
         rxn21890 rhea 15006
         rxn02157 1.14.11.13
         rxn21890 1.14.11.13
"""
rxn='rxn02157'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('RXN-115')
# if('rhea' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['rhea']=list()
# rxns_aliases_dict[rxn]['rhea'].append('15006')
if(rxn not in rxns_ecs_dict):
    rxns_ecs_dict[rxn]=list()
rxns_ecs_dict[rxn].append(' 1.14.11.13')

rxn='rxn21890'
# if('rhea' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['rhea']=list()
# rxns_aliases_dict[rxn]['rhea'].append('15006')
if(rxn not in rxns_ecs_dict):
    rxns_ecs_dict[rxn]=list()
rxns_ecs_dict[rxn].append(' 1.14.11.13')

del(rxns_aliases_dict['rxn42640'])
del(rxns_ecs_dict['rxn42640'])
"""
rxn46683
         rxn02158 MetaCyc RXN1F-170
         rxn02158 rhea 10105
         rxn26515 rhea 10105
         rxn02158 1.14.11.15
         rxn26515 1.14.11.15
"""
rxn='rxn02158'
# if('MetaCyc' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['MetaCyc']=list()
# rxns_aliases_dict[rxn]['MetaCyc'].append('RXN1F-170')
# if('rhea' not in rxns_aliases_dict[rxn]):
#    rxns_aliases_dict[rxn]['rhea']=list()
# rxns_aliases_dict[rxn]['rhea'].append('10105')
if(rxn not in rxns_ecs_dict):
    rxns_ecs_dict[rxn]=list()
rxns_ecs_dict[rxn].append(' 1.14.11.15')

rxn='rxn26515'
if('rhea' not in rxns_aliases_dict[rxn]):
    rxns_aliases_dict[rxn]['rhea']=list()
rxns_aliases_dict[rxn]['rhea'].append('10105')
if(rxn not in rxns_ecs_dict):
    rxns_ecs_dict[rxn]=list()
rxns_ecs_dict[rxn].append(' 1.14.11.15')

del(rxns_aliases_dict['rxn46683'])
del(rxns_ecs_dict['rxn46683'])

#reactions_helper.saveAliases(rxns_aliases_dict)
#reactions_helper.saveECs(rxns_ecs_dict)
