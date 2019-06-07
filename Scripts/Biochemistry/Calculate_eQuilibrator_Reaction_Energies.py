#!/usr/bin/env python
import logging
import json
import cobrakbase
from tqdm import tqdm
from BiochemPy.gibbs import ReactionEnergyCalculatorFromEquilibrator

modelseed_path = '/Users/fliu/workspace/jupyter/ModelSEEDDatabase/'

seed_id_to_equilibrator = {}

with open(modelseed_path + '/Biochemistry/Thermodynamics/seed_id_to_equilibrator.json', 'r') as f:
    seed_id_to_equilibrator = json.loads(f.read())

modelseed_local = cobrakbase.modelseed.from_local(modelseed_path)
    
rxn_calc_equilibrator = ReactionEnergyCalculatorFromEquilibrator(seed_id_to_equilibrator)

for seed_id in tqdm(modelseed_local.reactions):
    energy = None
    error = None
    metadata = {}
    status = 'OK'
    rxn = modelseed_local.get_seed_reaction(seed_id)
    if not rxn.is_obsolete or not rxn.is_transport:
        try:
            energy, error, metadata = rxn_calc_equilibrator.calculate_reaction_energy(rxn.cstoichiometry)
        except Exception as ex:
            status = ex
    else:
        status = 'obsolete/transport'
        
#TODO: format data and save it somewhere