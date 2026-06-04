#!/usr/bin/env python
import os,sys,json,glob
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions
from Estimate_Reaction_Reversibility import reversibility_from_energy

# dGPredictor (Wang et al. 2021, Metab Eng) predicts reaction dG directly from a
# group decomposition + ML model, output in kJ/mol. Predictions are staged as
# raw JSON in Biochemistry/Thermodynamics/dGPredictor/json_files/, keyed
# ModelSEED-rxn -> KEGG-R-id -> {dG_mean, dG_uncer}.
#
# This records the dGPredictor estimate ADDITIVELY: it stores a dGPredictor
# [energy, error, operator] record (kJ->kcal /4.184) in the JSON thermodynamics
# dict for EVERY reaction dGPredictor predicts, sitting next to the
# Group-Contribution / eQuilibrator records rather than replacing them. The
# operator is this estimate's own thermodynamic direction.
#
# It does NOT modify the canonical top-level deltag / deltagerr / reversibility:
# this is a description-only addition that leaves the served free-energy value
# under the control of the existing GC/eQ pipeline.

KJ_PER_KCAL = 4.184

label = "dGPredictor"
thermo_root = os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/dGPredictor/json_files/"

# rxn -> (kcal mean, kcal err) aggregated across its KEGG ids
dgp = dict()
for path in sorted(glob.glob(thermo_root+"reaction_*_dG.json")):
    with open(path) as fh:
        obj = json.load(fh)
    for rxn, kegg_map in obj.items():
        means=list(); uncs=list()
        if(isinstance(kegg_map, dict)):
            for kegg, payload in kegg_map.items():
                if(isinstance(payload, dict) and "dG_mean" in payload):
                    m = payload["dG_mean"]
                    u = payload.get("dG_uncer", 0.0)
                    if(isinstance(m,(int,float)) and m==m and abs(m)!=float("inf")):
                        means.append(m)
                        uncs.append(u if (isinstance(u,(int,float)) and u==u) else 0.0)
        if(means):
            dg_kcal = round(sum(means)/len(means)/KJ_PER_KCAL, 2)
            err_kcal = round(sum(uncs)/len(uncs)/KJ_PER_KCAL, 2)
            dgp[rxn] = (dg_kcal, err_kcal)

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

stored=0
for rxn in sorted(reactions_dict.keys()):
    robj = reactions_dict[rxn]

    if(robj.get('status') == 'EMPTY'):
        continue

    if(rxn not in dgp):
        continue

    (dg_kcal, err_kcal) = dgp[rxn]

    # ADDITIVE record: dGPredictor [energy, error, operator] alongside any
    # GC/eQ records. The operator is this estimate's own thermodynamic direction.
    operator = reversibility_from_energy(robj, dg_kcal, err_kcal)
    if(not isinstance(robj.get('thermodynamics'), dict)):
        robj['thermodynamics'] = dict()
    robj['thermodynamics'][label] = [dg_kcal, err_kcal, operator]
    stored+=1

print("dGPredictor reactions available: "+str(len(dgp)))
print("dGPredictor records stored additively: "+str(stored))
print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
