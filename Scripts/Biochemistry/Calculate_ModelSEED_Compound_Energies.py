#!/usr/bin/env python
import logging
import json
import cobrakbase
import pandas as pd
from BiochemPy.gibbs import get_delta_molanalysis

logger = logging.getLogger(__name__)

#read both 
modelseed_path = '/Users/fliu/workspace/jupyter/ModelSEEDDatabase/'
marvin_dg_file = '/Users/fliu/workspace/jupyter/python3/data/biochem/MolAnalysis_KEGG_MetaCyc_MarvinBeans_majorms_pH7.json'

kegg_compound_file = modelseed_path + '/Biochemistry/Aliases/Provenance/Primary_Databases/KEGG_Compounds.tbl'
meta_compound_file = modelseed_path + '/Biochemistry/Aliases/Provenance/Primary_Databases/MetaCyc_Compounds.tbl'


prov_kegg = pd.read_csv(kegg_compound_file, sep='\t', index_col=0)
prov_meta = pd.read_csv(meta_compound_file, sep='\t', index_col=0)

prov = {
    'KEGG' : prov_kegg,
    'MetaCyc' : prov_meta
}

modelseed_local = cobrakbase.modelseed.from_local(modelseed_path)

molanalysis_ph7 = None
with open(marvin_dg_file, 'r') as f:
    molanalysis_ph7 = json.loads(f.read())

databases = ['KEGG', 'MetaCyc']
seed_dg = {}

for seed_id in modelseed_local.compounds:
    (dG, dGe) = get_delta_molanalysis(seed_id, modelseed_local, molanalysis_ph7, prov, databases)
    if not None == dG and not None == dGe:
        seed_dg[seed_id] = (dG, dGe)
        logger.debug("[%s] %s %s", seed_id, dG, dGe)
    else:
        logger.debug("[%s] unable to get energy")

#TODO: format seed_dg and write it somewhere