#!/usr/bin/env python
import logging
import json
import cobrakbase
from tqdm import tqdm
import pandas as pd
from BiochemPy.gibbs import ReactionEnergyCalculatorFromCompound, CompoundEnergyCalculatorFromFile

modelseed_path = '/Users/fliu/workspace/jupyter/ModelSEEDDatabase/'

compound_energies = {}

#TODO: read the file from Thermodynamics/xxx
with open('/Users/fliu/workspace/jupyter/python3/data/biochem/gibbs/cpd_modelseed_ph7.json' , 'r') as f:
    compound_energies = json.loads(f.read())


modelseed_local = cobrakbase.modelseed.from_local(modelseed_path)

cpd_calc_from_file = CompoundEnergyCalculatorFromFile(compound_energies)
rxn_calc_modelseed = ReactionEnergyCalculatorFromCompound(cpd_calc_from_file)

for seed_id in tqdm(modelseed_local.reactions):
    energy = None
    error = None
    metadata = {}
    status = 'OK'
    rxn = modelseed_local.get_seed_reaction(seed_id)
    if not rxn.is_obsolete or not rxn.is_transport:
        try:
            energy, error, metadata = rxn_calc_modelseed.calculate_reaction_energy(rxn.cstoichiometry)
        except Exception as ex:
            status = ex
    else:
        status = 'obsolete/transport'
    #print(seed_id, energy, error, metadata, status)
    
#TODO: format data and save it somewhere