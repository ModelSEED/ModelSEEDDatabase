#!/usr/bin/env python
import os,sys,math
from equilibrator_api import ComponentContribution, Reaction, Q_
from BiochemPy import Compounds,Reactions

#We have to try and make sure that we use MetaNetX IDs for which an estimate of energy
#can be computed by eQuilibrator
thermodynamics_root=os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/"
file_name=thermodynamics_root+'eQuilibrator/MetaNetX_Compound_Energies.tbl'
Problem_Compounds=list()
with open(file_name) as file_handle:
    for line in file_handle.readlines():
        line = line.strip()
        array= line.split('\t')
        if('energy' not in array[1] and array[1] != 'nan'):
            continue
        Problem_Compounds.append(array[0])
file_handle.close()

#We have to map ModelSEED compound identifiers to the actual structures
#For which eQuilibrator generates an energy, so we compile the structures
#from eQuilibrator, skipping problematic ones
structures_root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
file_name=structures_root+'MetaNetX/Structures_in_ModelSEED_and_eQuilibrator.txt'
mnx_inchikey_dict=dict()
with open(file_name) as file_handle:
    for line in file_handle.readlines():
        line=line.strip()
        (mnx,inchikey)=line.split('\t')
        if(mnx in Problem_Compounds):
            continue

        if(inchikey not in mnx_inchikey_dict):
            mnx_inchikey_dict[inchikey]=mnx

        inchikey = "-".join(inchikey.split('-')[0:2])
        if(inchikey not in mnx_inchikey_dict):
            mnx_inchikey_dict[inchikey]=mnx

        inchikey = inchikey.split('-')[0]
        if(inchikey not in mnx_inchikey_dict):
            mnx_inchikey_dict[inchikey]=mnx

file_handle.close()

#Here we can cross-check the structures that are in ModelSEED to find ones where
#there is a match in eQuilibrator
compounds_helper = Compounds()
structures_dict = compounds_helper.loadStructures(["InChIKey"],["ModelSEED"])
seed_mnx_structural_map=dict()
for cpd in structures_dict:
    structure_type='InChIKey'
    if(structure_type not in structures_dict[cpd]):
        #The load structures function will return all compounds, so have to
        #Check that the structure is there
        continue

    #As these are unique structures, i.e. 1-1 mapping with compound id,
    #there's only ever one in each list for each compound
    structure = list(structures_dict[cpd][structure_type].keys())[0]

    #Here we check on three levels, we check the full string
    #Then the deprotonated string, then the structure alone
    #As per email from Elad and Moritz, we should not expect
    #Estimated energies to deviate between pseudoisomers (protons) and stereoisomers
    matched_mnx=None

    if(structure in mnx_inchikey_dict):
        matched_mnx=mnx_inchikey_dict[structure]
    else:
        structure = "-".join(structure.split('-')[0:2])
        if(structure in mnx_inchikey_dict):
            matched_mnx=mnx_inchikey_dict[structure]
        else:
            structure = structure.split('-')[0]
            if(structure in mnx_inchikey_dict):
                matched_mnx=mnx_inchikey_dict[structure]
    
    #As of 06/28/2019, there was:
    # 17,071 matches based on full inchikey
    # 17,863 matches using proton-neutral inchikey
    # 18,559 matches using stereo-neutral inchikey

    if(matched_mnx is None):
        continue

    seed_mnx_structural_map[cpd]=mnx_inchikey_dict[structure]

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

complete_mol_rxns_dict=dict()
incomplete_mol_rxns_dict=dict()
for rxn in reactions_dict:
    if(reactions_dict[rxn]['status']=='EMPTY'):
        continue

    rxn_cpds_array=reactions_helper.parseStoich(reactions_dict[rxn]["stoichiometry"])

    All_Mol=True
    Some_Mol=False
    for rgt in rxn_cpds_array:
        if(rgt['compound'] not in seed_mnx_structural_map):
            All_Mol=False
            Some_Mol=True

    if(All_Mol is True):
        complete_mol_rxns_dict[rxn]=1
    elif(Some_Mol is True):
        incomplete_mol_rxns_dict[rxn]=1

equilibrator_calculator = ComponentContribution(p_h=Q_(7.0), ionic_strength=Q_("0.25M"), temperature=Q_("298.15K"))
output_name=thermodynamics_root+'eQuilibrator/MetaNetX_Reaction_Energies.tbl'
output_handle=open(output_name,'w')
for rxn in reactions_dict:
    if(reactions_dict[rxn]['status']=="EMPTY"):
        continue

    notes_list=reactions_dict[rxn]['notes']
    if(not isinstance(notes_list,list)):
        notes_list=list()

    if(rxn not in complete_mol_rxns_dict):

        #'EQ' means equilibrator approach to calculating energies
        #'P' means partial, as in some of the reagents have energies calculated thus
        if(rxn in incomplete_mol_rxns_dict):
            if('EQC' in notes_list):
                notes_list.remove('EQC')
            if('EQP' not in notes_list):
                notes_list.append('EQP')

        if(len(notes_list)==0):
            reactions_dict[rxn]['notes']="null"
        else:
            reactions_dict[rxn]['notes']=notes_list
        continue

    #'EQ' means equilibrator approach to calculating energies
    #'C' means complete, as in all of the reagents have energies calculated thus
    if('EQP' in notes_list):
        notes_list.remove('EQP')
    if('EQC' not in notes_list):
        notes_list.append('EQC')

    rxn_cpds_array=reactions_helper.parseStoich(reactions_dict[rxn]["stoichiometry"])

    lhs=dict()
    rhs=dict()
    for rgt in rxn_cpds_array:
        if(rgt['compound'] not in seed_mnx_structural_map):
            OK = False
        else:
            mnx_id = 'mnx:'+seed_mnx_structural_map[rgt['compound']]

            if(rgt['coefficient'] < 0):
                lhs[mnx_id]=math.fabs(rgt['coefficient'])
            elif(rgt['coefficient'] > 0):
                rhs[mnx_id]=math.fabs(rgt['coefficient'])

    equation_str = ' + '.join([f'{value} {key}' for key, value in lhs.items()]) + \
        " = " + \
        ' + '.join([f'{value} {key}' for key, value in rhs.items()])

    equilibrator_reaction = Reaction.parse_formula(equation_str)

    try:
        dG0_prime, uncertainty = equilibrator_calculator.standard_dg_prime(equilibrator_reaction)
        dG0_prime = str(dG0_prime.to('kilocal / mole').magnitude)
        uncertainty = str(uncertainty.to('kilocal / mole').magnitude)

        ln_RI = equilibrator_calculator.ln_reversibility_index(equilibrator_reaction)
        if not type(ln_RI) == float:
            ln_RI = str(ln_RI.magnitude)

        output_handle.write("\t".join([rxn,dG0_prime,uncertainty,ln_RI])+"\n")
    except:
        output_handle.write("\t".join([rxn,"Unable to retrieve energy"])+"\n")

    #These are all true, as per earlier condition
    #print(OK)

print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
