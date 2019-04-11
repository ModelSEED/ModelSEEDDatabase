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

arguments = list(sys.argv)
#Pop script name
arguments = arguments[1:]
if(len(arguments) != 2 or 'cpd' not in arguments[0] or 'cpd' not in arguments[1]):
    print("Error: script must be initiated with the identifiers of the two compounds to merge")
    sys.exit()

From=arguments[1]
To=arguments[0]
arguments=sorted(arguments)

if(To != arguments[0]):
    print("Error: compound identifiers must be used in order")
    print("\tThe first compound identifier (in order) should be the one that is retained: "+arguments[0])
    sys.exit()

Dry_Run=0

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

if(compounds_dict[To]['is_obsolete']==1):
    print("Error: compound to be retained ("+To+")is marked as obsolete")
    sys.exit()

#Update compound links for compound to merge To
linked_cpds = list()
if(compounds_dict[To]['linked_compound'] is not None and compounds_dict[To]['linked_compound'] != "null"):
    linked_cpds=compounds_dict[To]['linked_compound'].split(";")
if(From not in linked_cpds):
    linked_cpds.append(From)
else:
    print("Warning: "+From+" already in 'linked_compound' list for "+To)
compounds_dict[To]['linked_compound']=";".join(sorted(linked_cpds))

#Update compound links for compound to merge From
linked_cpds = list()
if(compounds_dict[From]['linked_compound'] is not None and compounds_dict[From]['linked_compound'] != "null"):
    linked_cpds=compounds_dict[From]['linked_compound'].split(";")
if(To not in linked_cpds):
    linked_cpds.append(To)
else:
    print("Warning: "+To+" already in 'linked_compound' list for "+From)
compounds_dict[From]['linked_compound']=";".join(sorted(linked_cpds))

#Update obsolete flag for compound to merge From
compounds_dict[From]['is_obsolete']=1

if(Dry_Run==0):
    print("Saving merger of "+From+" to "+To)
    compounds_helper.saveCompounds(compounds_dict)
