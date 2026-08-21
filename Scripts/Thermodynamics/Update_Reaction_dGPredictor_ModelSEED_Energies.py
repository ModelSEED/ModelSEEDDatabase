#!/usr/bin/env python
import os,sys,json
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions
from Estimate_Reaction_Reversibility import reversibility_from_energy

# dGPredictor-ModelSEED is the dGPredictor group-contribution model (Wang et al.
# 2021, Metab Eng) RETRAINED on the ModelSEED compound structures rather than the
# original KEGG-only set. The retrain re-decomposes every ModelSEED compound that
# carries a complete structure into atom-centered fragments (radius 1 & 2),
# expanding the group vocabulary, and refits the same BayesianRidge model on the
# same 4,001 experimental measurements remapped into ModelSEED ID space. It
# predicts dG (kJ/mol, pH 7 / I 0.25 M / 298.15 K) for 31,924 reactions, including
# ~11,400 that the original KEGG-based model could not reach.
#
# This records the retrained estimate as its OWN additive method, "dGPredictor-
# ModelSEED", stored as a [energy, error, operator] record (kJ->kcal /4.184) in
# the JSON thermodynamics dict. It sits NEXT TO the existing "dGPredictor"
# (original KEGG-based) record and the Group-Contribution / eQuilibrator records
# rather than replacing any of them. The operator is this estimate's own
# thermodynamic direction.
#
# It does NOT modify the canonical top-level deltag / deltagerr / reversibility,
# nor the original "dGPredictor" record: this is a description-only addition that
# leaves the served free-energy value under the control of the existing GC/eQ
# pipeline.
#
# Predictions are staged as raw JSON in
# Biochemistry/Thermodynamics/dGPredictor/modelseed_retrained_dG.json, keyed
# ModelSEED-rxn -> {dG_mean, dG_uncer} in kJ/mol.

KJ_PER_KCAL = 4.184

label = "dGPredictor-ModelSEED"
staged = os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/dGPredictor/modelseed_retrained_dG.json"

# rxn -> (kcal mean, kcal err)
dgp = dict()
with open(staged) as fh:
    obj = json.load(fh)
for rxn, payload in obj.items():
    if(not isinstance(payload, dict) or "dG_mean" not in payload):
        continue
    m = payload["dG_mean"]
    u = payload.get("dG_uncer", 0.0)
    if(not (isinstance(m,(int,float)) and m==m and abs(m)!=float("inf"))):
        continue
    if(not (isinstance(u,(int,float)) and u==u and abs(u)!=float("inf"))):
        u = 0.0
    dg_kcal = round(m/KJ_PER_KCAL, 2)
    err_kcal = round(u/KJ_PER_KCAL, 2)
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

    # ADDITIVE record: dGPredictor-ModelSEED [energy, error, operator] alongside
    # any GC/eQ/original-dGPredictor records. The operator is this estimate's own
    # thermodynamic direction.
    operator = reversibility_from_energy(robj, dg_kcal, err_kcal)
    if(not isinstance(robj.get('thermodynamics'), dict)):
        robj['thermodynamics'] = dict()
    robj['thermodynamics'][label] = [dg_kcal, err_kcal, operator]
    stored+=1

print("dGPredictor-ModelSEED reactions available: "+str(len(dgp)))
print("dGPredictor-ModelSEED records stored additively: "+str(stored))
print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
