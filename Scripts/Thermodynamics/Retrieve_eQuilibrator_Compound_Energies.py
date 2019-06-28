#!/usr/bin/env python
import os
from equilibrator_api import ComponentContribution, Q_, Reaction

equilibrator_calculator = ComponentContribution(p_h=Q_(7.0), ionic_strength=Q_("0.25M"), temperature=Q_("298.15K"))

structures_root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
thermodynamics_root=os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/"
file_name=structures_root+'MetaNetX/Structures_in_ModelSEED_and_eQuilibrator.txt'
output_name=thermodynamics_root+'eQuilibrator/MetaNetX_Compound_Energies.tbl'
mnx_inchikey_dict=dict()
with open(file_name) as file_handle:
    for line in file_handle.readlines():
        line=line.strip()
        (mnx,inchikey)=line.split('\t')

        equilibrator_reaction = Reaction.parse_formula(' = ' + 'mnx:'+mnx)

        try:
            dG0_prime, uncertainty = equilibrator_calculator.standard_dg_prime(equilibrator_reaction)
            dG0_prime = str(dG0_prime.to('kilocal / mole').magnitude)
            uncertainty = str(uncertainty.to('kilocal / mole').magnitude)

            print("\t".join([mnx,dG0_prime,uncertainty]))
        except:
            print("\t".join([mnx,"Unable to retrieve energy"]))
