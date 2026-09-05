#!/usr/bin/env python

if __name__ == "__main__":
    # Validate arguments BEFORE importing anything or touching the database.
    # These scripts mutate the database, and without this an unknown flag or a
    # mistyped mode was silently ignored and the script ran with its defaults:
    # asking Estimate_Reaction_Reversibility.py for --help rewrote 122 files.
    # Placed above the imports so --help works even where a dependency is
    # missing from the path.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


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
# The operator is each estimate's OWN thermodynamic direction, computed with the
# rule set belonging to that method (eQuilibrator heuristics for eQuilibrator,
# Group-Contribution heuristics for everything else) applied to that method's dG.
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
        # Pass the label: each source is scored with its own rule set
        # (eQuilibrator gets the EQ heuristics, everything else GC).
        operator = reversibility_from_energy(rxn_obj, dg_val, dge_val, source=label)
        # Preserve anything past the operator. dGPredictor records carry a
        # fourth element (coverage); rebuilding the list as a fixed triple
        # would silently drop it on every backfill.
        new_values = [dg_val, dge_val, operator] + list(values[3:])
        entries += 1
        if(new_values != values):
            updated += 1
        thermo[label] = new_values

print("Per-method thermodynamics entries seen: " + str(entries))
print("Entries refreshed/added: " + str(updated))
print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
