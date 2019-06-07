#!/usr/bin/env python
import logging
import json
import cobrakbase
from tqdm import tqdm
from BiochemPy.gibbs import CompoundEnergyCalculatorFromEquilibrator

modelseed_path = '/Users/fliu/workspace/jupyter/ModelSEEDDatabase/'

seed_id_to_equilibrator = {}

with open(modelseed_path + '/Biochemistry/Thermodynamics/seed_id_to_equilibrator.json', 'r') as f:
    seed_id_to_equilibrator = json.loads(f.read())
    
modelseed_local = cobrakbase.modelseed.from_local(modelseed_path)

cpd_calc_equilibrator = energy.CompoundEnergyCalculatorFromEquilibrator(seed_id_to_equilibrator)

for seed_id in tqdm(modelseed_local.compounds, 'Compounds'):
    energy = None
    error = None
    metadata = {}
    status = 'OK'
    try:
        energy, error, metadata = cpd_calc_equilibrator.calculate_compound_energy(seed_id)
    except Exception as ex:
        status = ex
    #print(energy, error, metadata, status)
    
#TODO: format data and save it somewhere