#!/bin/bash
./Update_Compound_GroupContribution_Energies.py
./Update_Reaction_GroupContribution_Energies.py
./Estimate_Reaction_Reversibility.py GC
./Update_Compound_eQuilibrator_Energies.py
./Update_Reaction_eQuilibrator_Energies.py
./Estimate_Reaction_Reversibility.py EQ
# Record dGPredictor additively for every predicted reaction (no canonical change)
./Update_Reaction_dGPredictor_Energies.py
# Record the ModelSEED-retrained dGPredictor as its own additive method
./Update_Reaction_dGPredictor_ModelSEED_Energies.py
# Backfill/refresh the per-method [energy, error, operator] thermodynamics triples
./Add_Reaction_Thermodynamics_Operators.py
