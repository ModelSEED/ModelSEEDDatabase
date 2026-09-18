#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
cpds_aliases_dict = compounds_helper.loadMSAliases()
cpds_names_dict = compounds_helper.loadNames()

for cpd in compounds_dict:

    merged_aliases=dict()
    merged_names=list()

    if(cpd in cpds_aliases_dict):
        for source in cpds_aliases_dict[cpd]:
            if(source not in merged_aliases):
                merged_aliases[source]=list()
            for alias in cpds_aliases_dict[cpd][source]:
                if(alias not in merged_aliases[source]):
                    merged_aliases[source].append(alias)
    
    if(cpd in cpds_names_dict):
        for name in cpds_names_dict[cpd]:
            if(name not in merged_names):
                merged_names.append(name)
    
    cpds_aliases_dict[cpd]=merged_aliases
    cpds_names_dict[cpd]=merged_names
        
compounds_helper.saveNames(cpds_names_dict)
compounds_helper.saveAliases(cpds_aliases_dict)

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
rxns_aliases_dict = reactions_helper.loadMSAliases()
rxns_names_dict = reactions_helper.loadNames()
rxns_ecs_dict = reactions_helper.loadECs()

for rxn in reactions_dict:

    merged_aliases=dict()
    merged_names=list()
    merged_ecs=list()

    if(rxn in rxns_aliases_dict):
        for source in rxns_aliases_dict[rxn]:
            if(source not in merged_aliases):
                merged_aliases[source]=list()
            for alias in rxns_aliases_dict[rxn][source]:
                if(alias not in merged_aliases[source]):
                    merged_aliases[source].append(alias)
    
    if(rxn in rxns_names_dict):
        for name in rxns_names_dict[rxn]:
            if(name not in merged_names):
                merged_names.append(name)

    if(rxn in rxns_ecs_dict):
        for ec in rxns_ecs_dict[rxn]:
            if(ec not in merged_ecs):
                merged_ecs.append(ec)
    
    rxns_aliases_dict[rxn]=merged_aliases
    rxns_names_dict[rxn]=merged_names
    rxns_ecs_dict[rxn]=merged_ecs
        
reactions_helper.saveNames(rxns_names_dict)
reactions_helper.saveAliases(rxns_aliases_dict)
reactions_helper.saveECs(rxns_ecs_dict)

