#!/usr/bin/env python
"""Write reaction thermodynamics from the regenerated eQuilibrator energies.

Source: ``Biochemistry/Thermodynamics/eQuilibrator/ModelSEED_Reaction_Energies.tsv``,
computed directly from the structures ModelSEED holds, against an eQuilibrator
compound cache whose ``seed:`` accessions were repointed at those structures
(Path A) and whose training ``kegg:`` accessions were repointed so that
``kegg:X`` and ``seed:cpd#####`` resolve to the same compound (Path B), with
component-contribution retrained on the de-duplicated TECRDB
(``cc_params_dedup.npz``). Conditions are recorded in the file's header line.

This replaces the MetaNetX-mediated retrieval that fed
``MetaNetX_Reaction_Energies.tbl``. That path matched ModelSEED compounds to
MetaNetX ids down a progressively loosening InChIKey ladder -- full key, then
key minus protonation, then connectivity alone -- so a stored energy could rest
on a different stereoisomer or protonation state than the structure ModelSEED
actually holds. The .tbl is kept as the archived input of that pipeline (and is
still read by the reversibility-index tooling and its tests), but it no longer
feeds the database.

Coverage drops 25,028 -> 21,853, and most of that is deliberate. Of the ~6,200
reactions losing a value, ~4,000 were outside the component-contribution span
and ~1,900 do not balance -- the old code reported a number for both without
checking. Roughly 270 are a genuine cache gap. Reactions with no energy have
the key REMOVED rather than left holding a superseded number."""
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions
import _thermo_helpers as th

LABEL = 'eQuilibrator'

eq_reactions = th.parse_modelseed_energy_table(
    th.thermo_path('eQuilibrator', 'ModelSEED_Reaction_Energies.tsv'),
    id_col='reaction_id',
    dg_col='dg_prime_kcal_per_mol',
    err_col='uncertainty_kcal_per_mol')

print("%d reactions with an eQuilibrator energy" % len(eq_reactions))

th.run_reaction_table_update(Reactions(), LABEL, eq_reactions)
