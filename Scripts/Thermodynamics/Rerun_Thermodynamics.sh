#!/bin/bash
./Update_Compound_GroupContribution_Energies.py
./Update_Reaction_GroupContribution_Energies.py
./Estimate_Reaction_Reversibility.py GC
./Update_Compound_eQuilibrator_Energies.py
./Update_Reaction_eQuilibrator_Energies.py
./Estimate_Reaction_Reversibility.py EQ
