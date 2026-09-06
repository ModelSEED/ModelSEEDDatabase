#!/usr/bin/env python
import json
import os
import sys

sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions
from Estimate_Reaction_Reversibility import reversibility_from_energy

# dGPredictor (originally Wang et al. 2021, PLOS Comput Biol,
# doi:10.1371/journal.pcbi.1009448) — RETRAINED on the
# ModelSEED structure set with current-RDKit fragment canonicalization
# (2026-08-25). See Biochemistry/Thermodynamics/dGPredictor/README.md for
# training-pipeline details and the archived Wang-2020 baseline
# (dGPredictor-2020.tsv).
#
# Retrained predictions are staged at
# Biochemistry/Thermodynamics/dGPredictor/retrained_dG.json,
# schema: {rxn_id: {dG_mean_kJ_per_mol, dG_uncer_kJ_per_mol, coverage,
#                   frags_used, frags_oov}}
# where `coverage` is the fraction of the reaction's r=2 fragments present
# in the training vocabulary (0-1). Predictions with coverage below
# COVERAGE_FLOOR are not exported into the reactions JSON — they lack
# training-data support and would carry misleading uncertainty.
#
# This records the dGPredictor estimate ADDITIVELY: it stores a
# [energy, error, operator] triple (kJ->kcal /4.184) in the JSON
# thermodynamics dict for every reaction that passes the filter, next
# to (never replacing) the Group Contribution / eQuilibrator records.
# It does NOT modify canonical top-level deltag / deltagerr / reversibility.

KJ_PER_KCAL = 4.184
# No coverage floor. Measured against eQuilibrator on the 20,691 reactions both
# methods score, the reported uncertainty is well calibrated at EVERY coverage
# level -- the fraction within one combined sigma is 67-79% against 68% nominal,
# and never understates error. Coverage is also a poor proxy for accuracy:
# Spearman against |deviation| is -0.349 for coverage but +0.650 for sigma. The
# old 0.9 floor admitted 1,947 high-sigma predictions (median deviation 43.3
# kJ/mol) while rejecting 5,038 low-sigma ones (median deviation 12.3).
#
# Every prediction is therefore stored with its own uncertainty; consumers
# filter on sigma. Coverage is recorded as provenance, not as a filter.
COVERAGE_FLOOR = float(os.environ.get("DGP_COVERAGE_FLOOR", 0.0))

label = "dGPredictor"
staged = os.path.join(os.path.dirname(__file__),
                      "..", "..", "Biochemistry",
                      "Thermodynamics", "dGPredictor",
                      "retrained_dG.json")

with open(staged) as fh:
    payload = json.load(fh)

# rxn -> (kcal mean, kcal err); filter by coverage upfront.
dgp = dict()
n_below = 0
for rxn, entry in payload.items():
    if not isinstance(entry, dict) or "dG_mean" not in entry:
        continue
    if entry.get("coverage", 0) < COVERAGE_FLOOR:
        n_below += 1
        continue
    m = entry["dG_mean"]
    u = entry.get("dG_uncer", 0.0)
    if not (isinstance(m, (int, float)) and m == m and abs(m) != float("inf")):
        continue
    if not (isinstance(u, (int, float)) and u == u and abs(u) != float("inf")):
        u = 0.0
    dgp[rxn] = (round(m / KJ_PER_KCAL, 2), round(u / KJ_PER_KCAL, 2),
                round(entry.get("coverage", 0.0), 3))

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

stored = 0
for rxn in sorted(reactions_dict.keys()):
    robj = reactions_dict[rxn]
    if robj.get('status') == 'EMPTY':
        continue
    if rxn not in dgp:
        continue

    dg_kcal, err_kcal, cov = dgp[rxn]
    # source= picks the rule set; omitting it silently scores a
    # dGPredictor energy with Group-Contribution rules.
    operator = reversibility_from_energy(robj, dg_kcal, err_kcal,
                                         source=label)
    if not isinstance(robj.get('thermodynamics'), dict):
        robj['thermodynamics'] = dict()
    # [dg, err, operator, coverage] -- coverage appended, so positional
    # readers of [0]/[1]/[2] are unaffected.
    robj['thermodynamics'][label] = [dg_kcal, err_kcal, operator, cov]
    stored += 1

print(f"dGPredictor reactions available (coverage >= {COVERAGE_FLOOR}): {len(dgp)}")
print(f"dGPredictor predictions skipped for low coverage: {n_below}")
print(f"dGPredictor records stored additively: {stored}")
print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
