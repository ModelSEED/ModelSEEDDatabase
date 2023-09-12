#!/usr/bin/env python
import os,sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
structures_dict = compounds_helper.loadStructures(["SMILE","InChIKey"],["ModelSEED"])

############################################################################
##
## We apply/overwrite with eQuilibrator Energies
##
############################################################################

thermodynamics_root=os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/"
file_name=thermodynamics_root+'eQuilibrator/MetaNetX_Compound_Energies.tbl'
eq_compounds=dict()
with open(file_name) as file_handle:
    for line in file_handle.readlines():
        line = line.strip()
        array= line.split('\t')

        if('energy' in array[1] or array[1] == 'nan'):
            continue

        eq_compounds[array[0]]=[float("{0:.2f}".format(float(array[1]))),
                                float("{0:.2f}".format(float(array[2])))]
file_handle.close()

# print(len(eq_compounds))
# 19,432/22,391 (87%) MetaNetX records for which we can retrieve a compound formation energy

structures_root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
file_name=structures_root+'MetaNetX/Structures_in_ModelSEED_and_eQuilibrator.txt'
struct_mnx_dict=dict()
with open(file_name) as file_handle:
    for line in file_handle.readlines():
        line=line.strip()
        (mnx,inchikey)=line.split('\t')
    
        #For searching purposes we lose the protonation indicator
        inchikey="-".join(inchikey.split('-')[0:2])
        if(inchikey not in struct_mnx_dict):
            struct_mnx_dict[inchikey]=list()
        if(mnx not in struct_mnx_dict[inchikey]):
            struct_mnx_dict[inchikey].append(mnx)

file_handle.close()

# print(len(struct_mnx_dict))
# 18,206/19,432 (94%) MetaNetX records for which there is a unique structure

for cpd in sorted (compounds_dict.keys()):

    #Default energy and error
    dg_dge_list=[10000000.0,10000000.0]

    # Condition 1, no structure, use default
    # Condition 2, structure is InChIKey or SMILE

    structure_type=None
    if(cpd in structures_dict):
        if('InChIKey' in structures_dict[cpd]):
            structure_type = 'InChIKey' 
        elif('SMILE' in structures_dict[cpd]):
            structure_type = 'SMILE'

    structure = None
    if(structure_type is not None):
        structure = list(structures_dict[cpd][structure_type].keys())[0]

    if(structure is not None):

        deprotonated_struct="-".join(structure.split('-')[0:2])

        if(deprotonated_struct not in struct_mnx_dict):
            continue

        #In case where multiple energies because of distribution of bonds
        #Take lowest energy as most likely result of equilibrium
        #If the lowest energy is the default energy (i.e. 10000000)
        #We will still save it
        for mnx_id in struct_mnx_dict[deprotonated_struct]:
            if(mnx_id not in eq_compounds):
                continue
            
            if(eq_compounds[mnx_id][0] < dg_dge_list[0]):
                dg_dge_list=eq_compounds[mnx_id]
        
    #Here we indicate that we use the equilibrator value
    # values always saved as list of energy and error
    if(compounds_dict[cpd]['thermodynamics'] == "null"):
        compounds_dict[cpd]['thermodynamics'] = dict()
    if('eQuilibrator' not in compounds_dict[cpd]['thermodynamics']):
        compounds_dict[cpd]['thermodynamics']['eQuilibrator']=list()
    compounds_dict[cpd]['thermodynamics']['eQuilibrator']=dg_dge_list

print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)