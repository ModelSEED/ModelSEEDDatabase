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
# The script is idempotent: entries that already carry an operator (length 3)
# are left unchanged, and entries with no usable energy receive the '?' operator.
# It only recomputes from values already stored in the database, so it needs no
# upstream GC / eQuilibrator / dGPredictor inputs and can be re-run at any time.

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

        # already has an operator -> leave untouched (idempotent)
        if(len(values) >= 3):
            entries += 1
            continue

        (dg_val, dge_val) = (values[0], values[1])
        operator = reversibility_from_energy(rxn_obj, dg_val, dge_val)
        thermo[label] = [dg_val, dge_val, operator]
        entries += 1
        updated += 1

print("Per-method thermodynamics entries seen: " + str(entries))
print("Entries given an operator: " + str(updated))
print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
