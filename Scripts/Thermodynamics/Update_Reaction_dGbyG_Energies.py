#!/usr/bin/env python
"""Write reaction thermodynamics from the dGbyG energy table.

Source: ``Biochemistry/Thermodynamics/dGbyG/ModelSEED_Reaction_Energies.tsv``,
produced by ``Generate_dGbyG_Energies.py`` from the structures ModelSEED holds:
dGbyG (Fan et al. 2025, Cell Systems, doi:10.1016/j.cels.2025.101393), its
shipped 100-head ensemble, at pH 7.0, I 0.25 M, pMg 14, 298.15 K. The dGbyG
commit, model and conditions are recorded in the file's header line; method and
caveats are in ``Biochemistry/Thermodynamics/dGbyG/README.md``.

ADDITIVE: stores ``[dg, err, operator]`` in kcal/mol under
``thermodynamics['dGbyG']``, next to (never replacing) the other sources, and
leaves the canonical deltag / deltagerr / reversibility untouched. Only
``status == 'ok'`` rows carry an energy; a reaction without one has the key
REMOVED, so a missing key means dGbyG declined it (the table says why).

The operator is computed at write time with the rule set registered for
'dGbyG' in ``reversibility_heuristics.SOURCE_HEURISTIC_SET`` -- the reversibility
index at one sigma that eQuilibrator and dGPredictor use. An unregistered label
would silently fall back to the Group-Contribution rules, so this refuses to
run without the registration."""

if __name__ == "__main__":
    # Argument guard -- see "The argument guard" in Scripts/README.md.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions
import _thermo_helpers as th
from reversibility_heuristics import heuristics_for_source, DGB_HEURISTICS

LABEL = 'dGbyG'

if heuristics_for_source(LABEL) is not DGB_HEURISTICS:
    sys.exit("'dGbyG' is not registered in reversibility_heuristics.SOURCE_HEURISTIC_SET; "
             "its operators would be scored with Group-Contribution rules")

dgbyg_reactions = th.parse_modelseed_energy_table(
    th.thermo_path('dGbyG', 'ModelSEED_Reaction_Energies.tsv'),
    id_col='reaction_id',
    dg_col='dg_prime_kcal_per_mol',
    err_col='uncertainty_kcal_per_mol')

print("%d reactions with a dGbyG energy" % len(dgbyg_reactions))

th.run_reaction_table_update(Reactions(), LABEL, dgbyg_reactions)
