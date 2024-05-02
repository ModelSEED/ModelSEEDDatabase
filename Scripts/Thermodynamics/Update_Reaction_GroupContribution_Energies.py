#!/usr/bin/env python
import os,sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds, Reactions

label = 'Group contribution'

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

gc_cpds_dict=dict()
for cpd in compounds_dict:
    cpd_obj = compounds_dict[cpd]
    if(label in cpd_obj['thermodynamics'] and cpd_obj['thermodynamics'][label][0] != 10000000):
        gc_cpds_dict[cpd]=1

        if(cpd_obj['is_obsolete']):
            for link in cpd_obj['linked_compound'].split(';'):
                if(link in gc_cpds_dict):
                    continue

                link_obj = compounds_dict[link]
                if(label in link_obj['thermodynamics'] and link_obj['thermodynamics'][label][0] != 10000000):
                    gc_cpds_dict[link]=1

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

complete_gc_rxns_dict=dict()
incomplete_gc_rxns_dict=dict()
for rxn in reactions_dict:
    if(reactions_dict[rxn]['status']=='EMPTY'):
        continue

    rxn_cpds_array=reactions_dict[rxn]["stoichiometry"]

    all_gc=True
    some_gc=False
    for rgt in rxn_cpds_array:
        if(rgt['compound'] not in gc_cpds_dict):
            all_gc=False
            some_gc=True

    if(all_gc is True):
        complete_gc_rxns_dict[rxn]=1
    elif(some_gc is True):
        incomplete_gc_rxns_dict[rxn]=1

for rxn in reactions_dict:

    if(rxn not in complete_gc_rxns_dict and rxn not in incomplete_gc_rxns_dict):
        continue

    dg_dge_list = [10000000.0,10000000.0]

    if(rxn in complete_gc_rxns_dict):

        # build deltaG of reaction
        rxn_cpds_array=reactions_dict[rxn]["stoichiometry"]

        dg_sum=0.0
        dge_sum=0.0
        for rgt in rxn_cpds_array:
            if(rgt['compound'] not in gc_cpds_dict):
                print("Warning: wrong reaction: "+rxn)

            (dg,dge) = compounds_dict[rgt['compound']]['thermodynamics'][label]

            dg_sum += ( dg * rgt['coefficient'] )
            dge_sum+= ( dge * rgt['coefficient'] )**2

        dg_sum  =float("{0:.2f}".format(dg_sum))
        dge_sum =float("{0:.2f}".format(dge_sum**0.5))
        dg_dge_list = (dg_sum,dge_sum)

    # values always saved as list of energy and error
    if(not isinstance(reactions_dict[rxn]['thermodynamics'],dict)):
        reactions_dict[rxn]['thermodynamics'] = dict()
    if(label not in reactions_dict[rxn]['thermodynamics']):
        reactions_dict[rxn]['thermodynamics'][label]=list()
    reactions_dict[rxn]['thermodynamics'][label]=dg_dge_list

print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
