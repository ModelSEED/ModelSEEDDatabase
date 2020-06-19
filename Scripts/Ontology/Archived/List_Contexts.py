#!/usr/bin/env python
import os, sys
from BiochemPy import Compounds, InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

names=dict()
for cpd in compounds_dict:
    cpd_obj = compounds_dict[cpd]
    if(cpd_obj['aliases'] == "null"):
        continue

    if(cpd_obj['inchikey'] == ""):
        continue

    for entry in cpd_obj['aliases']:
        if('Name' not in entry):
            continue
        (type,names)=entry.split(': ',1)
        for name in names.split('; '):
            if("uino" in name.lower()):
                print(cpd_obj['id'],name)
