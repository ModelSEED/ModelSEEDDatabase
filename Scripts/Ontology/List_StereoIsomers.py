#!/usr/bin/env python
import os, sys
from BiochemPy import Compounds, InChIs

compounds_helper = Compounds()
structures_dict = compounds_helper.loadStructures(["InChI"],["ModelSEED"])
compounds_dict = compounds_helper.loadCompounds()

substructures=dict()
for cpd in compounds_dict:
    cpd_obj = compounds_dict[cpd]
    if(cpd_obj['is_obsolete']==0 and cpd_obj['inchikey'] != ""):
        substructure = cpd_obj['inchikey'].split('-')[0]
        if(substructure not in substructures):
            substructures[substructure]=list()
        substructures[substructure].append(cpd_obj['id'])

for substructure in sorted(substructures, key = lambda x: len(substructures[x]), reverse=True):
    print(substructure,len(substructures[substructure]), "|".join(substructures[substructure]))
