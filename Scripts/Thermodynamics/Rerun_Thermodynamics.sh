#!/bin/bash
./Update_Compound_GroupFormation_Energies.py
./Update_Reaction_GroupFormation_Energies.py
./Estimate_Reaction_Reversibility.py GF
./Update_Compound_eQuilibrator_Energies.py
./Update_Reaction_eQuilibrator_Energies.py
./Estimate_Reaction_Reversibility.py EQ
