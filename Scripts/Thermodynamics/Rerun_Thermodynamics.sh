#!/bin/bash
# Full thermodynamics regeneration, in dependency order.
#
# Every step reads a committed input under Biochemistry/Thermodynamics/ and is
# idempotent: run it twice against unchanged inputs and the second run produces
# no diff. That property is the check -- a non-empty diff on a second run means
# an input moved or a script is non-deterministic.
set -euo pipefail

# --- Group contribution -----------------------------------------------------
# Jankowski 2008 group energies via MFAToolkit, under Chris Henry's Convention A
# (H+ = -9.5 kcal/mol, H included in compound dGf). Inputs are the four
# MolAnalysis tables in Biochemistry/Thermodynamics/ModelSEED/.
./Update_Compound_GroupContribution_Energies.py
./Update_Reaction_GroupContribution_Energies.py

# --- eQuilibrator -----------------------------------------------------------
# Component contribution computed from ModelSEED's own structures. Inputs are
# ModelSEED_{Compound,Reaction}_Energies.tsv, which superseded the
# MetaNetX-mediated retrieval; see Biochemistry/Thermodynamics/eQuilibrator/README.md.
./Update_Compound_eQuilibrator_Energies.py
./Update_Reaction_eQuilibrator_Energies.py

# --- dGPredictor ------------------------------------------------------------
# Retrained on the de-duplicated TECRDB, the same measurements
# component-contribution is fitted to. Every prediction is installed with its
# own uncertainty; there is no coverage floor. See
# Biochemistry/Thermodynamics/dGPredictor/README.md.
./Update_Reaction_dGPredictor_Energies.py
./Update_Compound_dGPredictor_Energies.py

# --- Per-source direction operators -----------------------------------------
# Each updater already computes its own operator at write time, using that
# source's rule set. This is a consistency backfill: on a clean run it reports
# "Entries refreshed/added: 0". A non-zero count means some record was written
# with the wrong rule set.
./Add_Reaction_Thermodynamics_Operators.py

# --- Canonical fields (deltag / deltagerr / reversibility) ------------------
# NOT run above, deliberately. These write the top-level canonical fields,
# which are slated for removal in favour of the additive per-source dict. They
# are listed here so the full pipeline stays documented:
#
#   ./Estimate_Reaction_Reversibility.py GC
#   ./Estimate_Reaction_Reversibility.py EQ
#   ./Promote_Reaction_Thermodynamics_to_Canonical.py
#
# Promotion must run LAST if it runs at all -- it weighs reversibility when
# choosing which source to promote, so running it between source updates
# promotes a half-regenerated picture. Note also that it never overwrites an
# existing canonical deltag, so on a database that already has one it is close
# to a no-op.
