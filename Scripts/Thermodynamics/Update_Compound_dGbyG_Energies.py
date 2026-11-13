#!/usr/bin/env python
"""Write compound formation energies from the dGbyG energy table.

Source: ``Biochemistry/Thermodynamics/dGbyG/ModelSEED_Compound_Energies.tsv``
-- same run, same model and same conditions as the reaction table; see
``Update_Reaction_dGbyG_Energies.py`` and the README beside the table for the
provenance.

ADDITIVE: stores ``[dg, err]`` in kcal/mol under ``thermodynamics['dGbyG']``
next to the other sources. Only ``status == 'ok'`` rows carry an energy;
compounds without a structure, with an R group ('*'), or that dGbyG could not
parse have the key REMOVED.

dGbyG's formation energies sit on dGbyG's own reference. Against eQuilibrator
(the 16,651 live compounds where both carry a usable value) they correlate
closely, Pearson r 0.995, but differ compound by compound by a median 8.0
kcal/mol -- an offset that hydrogen count and charge do not explain (R^2 0.002)
and that largely cancels in balanced reactions. Compare sources at the reaction
level.

Compounds store ``[dg, dge]``: a formation energy has no direction."""

if __name__ == "__main__":
    # Argument guard -- see "The argument guard" in Scripts/README.md.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds
import _thermo_helpers as th

LABEL = 'dGbyG'

dgbyg_compounds = th.parse_modelseed_energy_table(
    th.thermo_path('dGbyG', 'ModelSEED_Compound_Energies.tsv'),
    id_col='compound_id',
    dg_col='dgf_prime_kcal_per_mol',
    err_col='uncertainty_kcal_per_mol')

print("%d compounds with a dGbyG formation energy" % len(dgbyg_compounds))

th.run_compound_table_update(Compounds(), LABEL, dgbyg_compounds)
