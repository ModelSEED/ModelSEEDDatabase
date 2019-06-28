#!/usr/bin/env python
import os,sys
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

        eq_compounds[array[0]]={'dg':"{0:.2f}".format(float(array[1])),'dge':"{0:.2f}".format(float(array[2]))}

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
        
        #This check makes sure that we use MetaNetX IDs for which we know we
        #can retrieve energies
        if(mnx not in eq_compounds):
            continue

        if('struct' not in eq_compounds[mnx]):
            eq_compounds[mnx]['struct']=inchikey

            #For searching purposes we lose the protonation indicator
            inchikey="-".join(inchikey.split('-')[0:2])
            struct_mnx_dict[inchikey]=mnx

file_handle.close()

# print(len(struct_mnx_dict))
# 18,206/19,432 (94%) MetaNetX records for which there is a unique structure

seed_mnx_map=dict()
for cpd in structures_dict:
    structure_type='InChIKey'
    if(structure_type not in structures_dict[cpd]):
        continue

    structure = list(structures_dict[cpd][structure_type].keys())[0]
    dp_struct="-".join(structure.split('-')[0:2])

    if(dp_struct not in struct_mnx_dict):
        continue

    seed_mnx_map[cpd]=eq_compounds[struct_mnx_dict[dp_struct]]

# print(len(seed_mnx_map))
# 17,863 ModelSEED compounds assigned eQuilibrator energies

file_handle = open('GroupFormation_eQuilibrator_Comparison.txt', 'w')
file_handle.write('ID\tGF\tEQ\n')
for cpd in sorted (compounds_dict.keys()):

    cpd_gf='nan'
    cpd_eq='nan'

    if("GF" in compounds_dict[cpd]['notes']):
        cpd_gf='|'.join([str(compounds_dict[cpd]['deltag']),str(compounds_dict[cpd]['deltagerr'])])

    if(cpd in seed_mnx_map):
        cpd_eq='|'.join([seed_mnx_map[cpd]['dg'],seed_mnx_map[cpd]['dge']])

    #Write the values to file. I'm not making exception, but I'm
    #including this `if` statement as a comment to remind me of ways in which
    #they can't be directly compared
    #if(compounds_dict[cpd]['deltag']!=10000000 and cpd_gf != 'nan' and cpd_eq != 'nan'):
    file_handle.write('\t'.join([cpd,cpd_gf,cpd_eq])+'\n')

    #Having printed to file, we skip where there is no eQuilibrator estimate
    if(cpd not in seed_mnx_map):
        continue

    #Here we establish an arbitrary threshold of 50 for the error, if the error
    #is too big, we don't use it

    if(float(seed_mnx_map[cpd]['dge']) > 50):
        continue

    if(compounds_dict[cpd]['notes'] == "null"):
        compounds_dict[cpd]['notes']=list()

    compounds_dict[cpd]['deltag']=float(seed_mnx_map[cpd]['dg'])
    compounds_dict[cpd]['deltagerr']=float(seed_mnx_map[cpd]['dge'])
    if('EQ' not in compounds_dict[cpd]['notes']):
        compounds_dict[cpd]['notes'].append('EQ')

file_handle.close()
print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
