#!/usr/bin/env python
import os
from equilibrator_api import ComponentContribution, Q_, Reaction, ccache

equilibrator_calculator = ComponentContribution(p_h=Q_(7.0), ionic_strength=Q_("0.25M"), temperature=Q_("298.15K"))

structures_root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
thermodynamics_root=os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/"
file_name=structures_root+'MetaNetX/Structures_in_ModelSEED_and_eQuilibrator.txt'
output_name=thermodynamics_root+'eQuilibrator/MetaNetX_Compound_Energies.tbl'
output_handle=open(output_name,'w')
mnx_inchikey_dict=dict()
with open(file_name) as file_handle:
    for line in file_handle.readlines():
        line=line.strip()
        (mnx,inchikey)=line.split('\t')

        equilibrator_reaction = Reaction.parse_formula(ccache.get_compound, ' = ' + mnx)

        try:
            result = equilibrator_calculator.standard_dg_prime(equilibrator_reaction)
            dG0_prime = str(result.value.to('kilocal / mole').magnitude)
            uncertainty = str(result.error.to('kilocal / mole').magnitude)

            output_handle.write("\t".join([mnx,dG0_prime,uncertainty])+"\n")
        except Exception as e:
            output_handle.write("\t".join([mnx,"Unable to retrieve energy"])+"\n")
            print(e)
