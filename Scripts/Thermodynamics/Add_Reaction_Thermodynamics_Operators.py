#!/usr/bin/env python
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions
from Estimate_Reaction_Reversibility import reversibility_from_energy

# Backfill the per-method thermodynamic-direction operator onto every entry in a
# reaction's `thermodynamics` dict, turning the legacy [energy, error] pairs into
# the additive [energy, error, operator] triples, e.g.
#
#   "thermodynamics": {
#       "Group contribution": [4.15, 1.22, "="],
#       "eQuilibrator":       [-3.46, 0.05, ">"]
#   }
#
# The operator is each estimate's OWN thermodynamic direction (computed with the
# same heuristic as the canonical reversibility, applied to that method's dG).
# The canonical top-level deltag / deltagerr / reversibility fields are NEVER
# touched: these per-method records are added next to, not in place of, the
# existing values.
#
# The operator is (re)computed for every entry from its stored [energy, error]:
# legacy [energy, error] pairs are promoted to triples, and entries that already
# carry an operator are REFRESHED so they track the current cascade (run this
# after any change to Estimate_Reaction_Reversibility.py's heuristics). Entries
# with no usable energy receive the '?' operator. The recompute is deterministic,
# so the script is idempotent and needs no upstream GC / eQuilibrator /
# dGPredictor inputs -- it can be re-run at any time.

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

updated = 0
entries = 0
for rxn in sorted(reactions_dict.keys()):
    rxn_obj = reactions_dict[rxn]

    thermo = rxn_obj.get('thermodynamics')
    if(not isinstance(thermo, dict)):
        continue

    for label in thermo:
        values = thermo[label]
        if(not isinstance(values, list) or len(values) < 2):
            continue

        (dg_val, dge_val) = (values[0], values[1])
        operator = reversibility_from_energy(rxn_obj, dg_val, dge_val)
        new_values = [dg_val, dge_val, operator]
        entries += 1
        if(new_values != values):
            updated += 1
        thermo[label] = new_values

print("Per-method thermodynamics entries seen: " + str(entries))
print("Entries refreshed/added: " + str(updated))
print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
