#!/usr/bin/env python
import json
import os
import sys

sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds

# Per-compound ΔfG'° from the RETRAINED dGPredictor model (2026-08-25).
# See Biochemistry/Thermodynamics/dGPredictor/README.md for training-pipeline
# details.
#
# The retrained model was fit on the TECRDB training set alongside 224
# formation-energy pseudo-reactions of the form `∅ → cpd`, so it has learned
# per-compound ΔfG'° directly. Predictions are staged at
# Biochemistry/Thermodynamics/dGPredictor/retrained_dG_compounds.json,
# schema: {cpd_id: {dG_mean_kJ_per_mol, dG_uncer_kJ_per_mol, coverage,
#                   frags_used, frags_oov}}
# where `coverage` is the fraction of the compound's r=2 fragments present
# in the training vocabulary (0-1).
#
# This records the dGPredictor estimate ADDITIVELY: it stores a
# [energy, error] pair (kJ->kcal /4.184) in the JSON thermodynamics dict for
# every compound that passes the filter, next to (never replacing) the
# Group Contribution / eQuilibrator records. It does NOT modify canonical
# top-level deltag / deltagerr.
#
# Compounds previously carried only 'Group contribution' and 'eQuilibrator';
# adding dGPredictor gives a third source for method comparison in the paper
# (Seaver et al. 2026 NAR Database Issue).

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
                      "retrained_dG_compounds.json")

with open(staged) as fh:
    payload = json.load(fh)

dgp = dict()
n_below = 0
for cpd, entry in payload.items():
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
    dgp[cpd] = (round(m / KJ_PER_KCAL, 2), round(u / KJ_PER_KCAL, 2),
                round(entry.get("coverage", 0.0), 3))

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

stored = 0
for cpd in sorted(compounds_dict.keys()):
    cobj = compounds_dict[cpd]
    if cpd not in dgp:
        continue
    dg_kcal, err_kcal, cov = dgp[cpd]
    if not isinstance(cobj.get('thermodynamics'), dict):
        cobj['thermodynamics'] = dict()
    # [dg, err, coverage] -- compounds carry no operator, a formation
    # energy having no direction.
    cobj['thermodynamics'][label] = [dg_kcal, err_kcal, cov]
    stored += 1

print(f"dGPredictor compounds available (coverage >= {COVERAGE_FLOOR}): {len(dgp)}")
print(f"dGPredictor predictions skipped for low coverage: {n_below}")
print(f"dGPredictor records stored additively: {stored}")
print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
