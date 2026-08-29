#!/usr/bin/env python
"""Write compound formation energies from the regenerated eQuilibrator run.

Source: ``Biochemistry/Thermodynamics/eQuilibrator/ModelSEED_Compound_Energies.tsv``
-- same cache, same parameters and same conditions as the reaction table; see
``Update_Reaction_eQuilibrator_Energies.py`` for the provenance.

The lookup is now by ModelSEED compound id, because the cache resolves
``seed:cpd#####`` to the structure ModelSEED holds for it. The previous script
had to work backwards: deprotonate the stored InChIKey to two segments, look
that up in a MetaNetX index, and take the lowest energy among the MetaNetX
records that matched. Every step of that was a guess about which molecule the
accession meant, which is the thing this regeneration exists to remove.

Coverage drops 30,607 -> 16,372. The compound side loses proportionally more
than the reaction side because 43.7% of ModelSEED compounds fall outside the
component-contribution span and another 20.5% are absent from the cache -- both
previously papered over by the loosening-InChIKey fallback. Compounds with no
energy have the key REMOVED rather than left holding a superseded number.

Compounds store ``[dg, dge]``: a formation energy has no direction."""
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds
import _thermo_helpers as th

LABEL = 'eQuilibrator'

eq_compounds = th.parse_modelseed_energy_table(
    th.thermo_path('eQuilibrator', 'ModelSEED_Compound_Energies.tsv'),
    id_col='compound_id',
    dg_col='dgf_prime_kcal_per_mol',
    err_col='uncertainty_kcal_per_mol')

print("%d compounds with an eQuilibrator formation energy" % len(eq_compounds))

th.run_compound_table_update(Compounds(), LABEL, eq_compounds)
