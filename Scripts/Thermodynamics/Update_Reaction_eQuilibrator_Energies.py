#!/usr/bin/env python
"""Overwrite reaction thermodynamics with eQuilibrator energies.

Source: ``Biochemistry/Thermodynamics/eQuilibrator/MetaNetX_Reaction_Energies.tbl``
(from ``Retrieve_eQuilibrator_Reactions_Energies.py``).

Reactions absent from the table are left untouched. There are reactions
for which we have eQuilibrator structures for every reagent (labelled
``EQC`` in reaction notes) yet still no estimated reaction energy —
likely because at least one compound's energy lookup failed."""
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions
import _thermo_helpers as th

LABEL = 'eQuilibrator'

eq_reactions = th.parse_two_col_energy_table(
    th.thermo_path('eQuilibrator', 'MetaNetX_Reaction_Energies.tbl'))

# ~13,874 ModelSEED reactions covered.

th.run_reaction_lookup_update(Reactions(), LABEL, eq_reactions)
