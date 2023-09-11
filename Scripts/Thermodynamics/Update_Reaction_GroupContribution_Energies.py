#!/usr/bin/env python
import os,sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds, Reactions

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

mol_cpds_dict=dict()
for cpd in compounds_dict:
    if('GC' in compounds_dict[cpd]['notes'] and compounds_dict[cpd]['deltag'] != 10000000):
        mol_cpds_dict[cpd]=1

        if(compounds_dict[cpd]['is_obsolete']):
            for link in compounds_dict[cpd]['linked_compound'].split(';'):
                if('GC' in compounds_dict[link]['notes'] and compounds_dict[cpd]['deltag'] != 10000000):
                    mol_cpds_dict[link]=1

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

complete_mol_rxns_dict=dict()
incomplete_mol_rxns_dict=dict()
for rxn in reactions_dict:
    if(reactions_dict[rxn]['status']=='EMPTY'):
        continue

    rxn_cpds_array=reactions_dict[rxn]["stoichiometry"]

    All_Mol=True
    Some_Mol=False
    for rgt in rxn_cpds_array:
        if(rgt['compound'] not in mol_cpds_dict):
            All_Mol=False
            Some_Mol=True

    if(All_Mol is True):
        complete_mol_rxns_dict[rxn]=1
    elif(Some_Mol is True):
        incomplete_mol_rxns_dict[rxn]=1

for rxn in reactions_dict:

    notes_list=reactions_dict[rxn]['notes']
    if(not isinstance(notes_list,list)):
        notes_list=list()

    if(rxn not in complete_mol_rxns_dict):

        #'GC' means group contribution approach to calculating energies
        #'P' means partial, as in some of the reagents have energies calculated thus
        if(rxn in incomplete_mol_rxns_dict):
            if('GCC' in notes_list):
                notes_list.remove('GCC')
            if('GCP' not in notes_list):
                notes_list.append('GCP')

        reactions_dict[rxn]['deltag']=10000000.0
        reactions_dict[rxn]['deltagerr']=10000000.0

        if(len(notes_list)==0):
            reactions_dict[rxn]['notes']="null"
        else:
            reactions_dict[rxn]['notes']=notes_list
        continue

    #'GC' means group contribution approach to calculating energies
    #'C' means complete, as in all of the reagents have energies calculated thus
    if('GCP' in notes_list):
        notes_list.remove('GCP')
    if('GCC' not in notes_list):
        notes_list.append('GCC')

    rxn_cpds_array=reactions_dict[rxn]["stoichiometry"]

    #thermodynamics
    dg_sum=0.0
    dge_sum=0.0
    for rgt in rxn_cpds_array:
        if(rgt['compound'] not in mol_cpds_dict):
            print("Warning: wrong reaction: "+rxn)

        dg_sum+= ( compounds_dict[rgt['compound']]['deltag'] * rgt['coefficient'] )
        dge_sum+= ( compounds_dict[rgt['compound']]['deltagerr'] * rgt['coefficient'] )**2

    dg_sum="{0:.2f}".format(dg_sum)
    dge_sum = "{0:.2f}".format(dge_sum**0.5)

    reactions_dict[rxn]['deltag']=float(dg_sum)
    reactions_dict[rxn]['deltagerr']=float(dge_sum)
    reactions_dict[rxn]['notes']=notes_list

print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
