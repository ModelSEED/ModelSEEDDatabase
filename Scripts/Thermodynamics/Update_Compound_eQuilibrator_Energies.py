#!/usr/bin/env python
"""Overwrite compound thermodynamics with eQuilibrator energies.

Source: ``Biochemistry/Thermodynamics/eQuilibrator/MetaNetX_Compound_Energies.tbl``
(from ``Retrieve_eQuilibrator_Compound_Energies.py``) plus the MetaNetX
structure-to-InChIKey index under ``Biochemistry/Structures/MetaNetX/``.

Per-compound resolution:
  - structure picked via ``InChIKey`` then ``SMILE`` preference
  - structure deprotonated to a 2-segment InChIKey, then looked up against
    the MetaNetX index to recover candidate MNX ids
  - lowest dg among matching MNX records wins

Skips writing entirely when the compound has a structure that doesn't
appear in the MetaNetX index (preserves original eQuilibrator skip
semantics). Writes the default sentinel when no structure is available."""
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds
import _thermo_helpers as th

LABEL = 'eQuilibrator'

compounds_helper = Compounds()
eq_compounds = th.parse_two_col_energy_table(
    th.thermo_path('eQuilibrator', 'MetaNetX_Compound_Energies.tbl'))
struct_mnx_dict = th.parse_mnx_inchikey_index(
    th.structures_path('MetaNetX', 'Structures_in_ModelSEED_and_eQuilibrator.txt'))

# 19,432/22,391 (87%) MetaNetX records carry a compound formation energy.
# 18,206/19,432 (94%) of those map to a unique structure.


def resolve(cpd, stype, structure, aliases):
    deprotonated = "-".join(structure.split('-')[0:2])
    if deprotonated not in struct_mnx_dict:
        return None
    return th.lowest_energy_eq_style(struct_mnx_dict[deprotonated], eq_compounds)


th.run_compound_update(compounds_helper, LABEL, resolve, on_no_structure='default')
