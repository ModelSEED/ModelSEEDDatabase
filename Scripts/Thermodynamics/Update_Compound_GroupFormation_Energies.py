#!/usr/bin/env python
import os,sys
from BiochemPy import Compounds

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
structures_dict = compounds_helper.loadStructures(["SMILE","InChIKey"],["ModelSEED"])

############################################################################
##
## First we apply Group Formation Energies
##
############################################################################
#In the case where there was originally conflicting structures, we only want
#The energy for the structure that was curated, and we're using the provenance here
#To delete the aliases from structures_dict that are from conflicting structures
structures_root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
file_name=structures_root+'All_ModelSEED_Structures.txt'
all_structures_dict=dict()
with open(file_name) as file_handle:
    for line in file_handle.readlines():
        line=line.strip()
        array=line.split('\t')
        if(array[0] not in all_structures_dict):
            all_structures_dict[array[0]]=dict()
        if(array[7] not in all_structures_dict[array[0]]):
            all_structures_dict[array[0]][array[7]]=dict()
        all_structures_dict[array[0]][array[7]][array[3]]=1

for cpd in structures_dict:
    structure_type='InChIKey'
    if(structure_type not in structures_dict[cpd]):
        structure_type='SMILE'

    structure = list(structures_dict[cpd][structure_type].keys())[0]
    aliases_list = list(structures_dict[cpd][structure_type][structure]['alias'])
    new_aliases_list = list()
    for alias in aliases_list:
        if(alias in all_structures_dict[cpd][structure]):
            new_aliases_list.append(alias)
    structures_dict[cpd][structure_type][structure]['alias']=new_aliases_list
############################################################################

thermodynamics_root=os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/"
thermodynamics_dict=dict()
for source in ["KEGG","MetaCyc"]:
    for process in ["Charged","Original"]:
        file_name=thermodynamics_root+'ModelSEED/'+source+'_'+process+'_MolAnalysis.tbl'
        with open(file_name) as file_handle:
            for line in file_handle.readlines():
                line=line.strip()
                array=line.split('\t')
                if(array[0] not in thermodynamics_dict):
                    thermodynamics_dict[array[0]]={'dg':"{0:.2f}".format(float(array[7])),'dge':"{0:.2f}".format(float(array[8]))}
                else:
                    #There's a few (~20) cases where the protonated mol file had a 'NoGroup' cue added by MFAToolkit
                    #So using Original energy
                    if(thermodynamics_dict[array[0]]['dg'] == "10000000.00" and array[7] != "1e+07"):
                        thermodynamics_dict[array[0]]={'dg':"{0:.2f}".format(float(array[7])),'dge':"{0:.2f}".format(float(array[8]))}

for cpd in sorted (compounds_dict.keys()):

    #Energies computed from structures, if no structure, don't even _have_ energies
    if(cpd not in structures_dict):
        compounds_dict[cpd]['deltag']=10000000.0
        compounds_dict[cpd]['deltagerr']=10000000.0
        compounds_dict[cpd]['notes']="null"

    else:

        #If inchikey, use that, else use SMILE
        structure_type='InChIKey'
        if(structure_type not in structures_dict[cpd]):
            structure_type='SMILE'

        structure = list(structures_dict[cpd][structure_type].keys())[0]
        energies_dict=dict()
        for alias in structures_dict[cpd][structure_type][structure]['alias']:
            if(alias not in thermodynamics_dict):
                continue
            energies_dict[float(thermodynamics_dict[alias]['dg'])]=float(thermodynamics_dict[alias]['dge'])

        #In case where multiple energies because of distribution of bonds
        #Take lowest energy as most likely result of equilibrium
        lowest_dg=10000000.0
        lowest_dge=10000000.0
        for energy in energies_dict:
            if(energy < lowest_dg):
                lowest_dg=energy
                lowest_dge=energies_dict[energy]

        #If the lowest energy is the default energy (i.e. 10000000)
        #We will still save it, but we're taking notes

        compounds_dict[cpd]['deltag']=lowest_dg
        compounds_dict[cpd]['deltagerr']=lowest_dge
        compounds_dict[cpd]['notes']=['GF'] #Meaning group formation approach

############################################################################
##
## Second we apply/overwrite with eQuilibrator Energies
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

    compounds_dict[cpd]['deltag']=seed_mnx_map[cpd]['dg']
    compounds_dict[cpd]['deltagerr']=seed_mnx_map[cpd]['dge']
    compounds_dict[cpd]['notes'].append('EQ')

file_handle.close()
print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
