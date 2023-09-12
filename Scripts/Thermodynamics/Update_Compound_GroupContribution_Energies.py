#!/usr/bin/env python
import os,sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
structures_dict = compounds_helper.loadStructures(["SMILE","InChIKey"],["ModelSEED"])

############################################################################
##
## First we apply Group Contribution Energies
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

    #Default energy and error
    lowest_dg=10000000.0
    lowest_dge=10000000.0

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
        energies_dict=dict()
        for alias in structures_dict[cpd][structure_type][structure]['alias']:
            if(alias not in thermodynamics_dict):
                continue
            energies_dict[float(thermodynamics_dict[alias]['dg'])]=float(thermodynamics_dict[alias]['dge'])

        #In case where multiple energies because of distribution of bonds
        #Take lowest energy as most likely result of equilibrium
        #If the lowest energy is the default energy (i.e. 10000000)
        #We will still save it
        for energy in energies_dict:
            if(energy < lowest_dg):
                lowest_dg=energy
                lowest_dge=energies_dict[energy]

    # values always saved as list of energy and error
    if(compounds_dict[cpd]['thermodynamics'] == "null"):
        compounds_dict[cpd]['thermodynamics'] = dict()
    if('Group contribution' not in compounds_dict[cpd]['thermodynamics']):
        compounds_dict[cpd]['thermodynamics']['Group contribution']=list()
    compounds_dict[cpd]['thermodynamics']['Group contribution'].append(lowest_dg)
    compounds_dict[cpd]['thermodynamics']['Group contribution'].append(lowest_dge)

print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
