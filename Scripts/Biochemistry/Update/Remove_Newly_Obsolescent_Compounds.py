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

# We actually don't want obsolete reactions and compounds in our database
# So we're striving to remove any 'new' ones that are obsolete
# Any information attached to them should be associated with their linked counterpart
# We need to retain older compounds that are now obsolete as these may be present in prior published models

# The number used here is the last compound entered before we re-integrated updates from KEGG and MetaCyc
# In the fall of 2018, so after this point, we'll take out obsolete compounds
last_cpd_str='cpd31000'
last_cpd_int=int(last_cpd_str[3:])

delete_cpds=list()
for cpd in compounds_dict:
    cpd_int = int(cpd[3:])
    if(cpd_int > last_cpd_int and compounds_dict[cpd]['is_obsolete']):
        delete_cpds.append(cpd)

for cpd in delete_cpds:

    # I need to update the linked_compound field to make sure that the removed compound is deleted
    for lnkd_cpd in compounds_dict[cpd]['linked_compound'].split(';'):
        linked_compound_list = compounds_dict[lnkd_cpd]['linked_compound'].split(';')
        linked_compound_list.remove(cpd)
        compounds_dict[lnkd_cpd]['linked_compound']=linked_compound_list

        # I need to move the names, aliases, and ec numbers to linked compounds
        # This should already have been done when merging compounds, but doing it here to double-check
        if(cpd in cpds_aliases_dict):
            for source in cpds_aliases_dict[cpd]:
                if(source not in cpds_aliases_dict[lnkd_cpd]):
                    cpds_aliases_dict[lnkd_cpd][source]=list()
                for alias in cpds_aliases_dict[cpd][source]:
                    if(alias not in cpds_aliases_dict[lnkd_cpd][source]):
                        cpds_aliases_dict[lnkd_cpd][source].append(alias)
                        print("Warning: adding "+alias+" for "+source+" to "+lnkd_cpd)

        if(cpd in cpds_names_dict):
            for name in cpds_names_dict[cpd]:
                if(name not in cpds_names_dict[lnkd_cpd]):
                    cpds_names_dict[lnkd_cpd].append(name)
                    print("Warning: adding "+name+" to "+lnkd_cpd)

    del(compounds_dict[cpd])
    if(cpd in cpds_aliases_dict):
        del(cpds_aliases_dict[cpd])
    if(cpd in cpds_names_dict):
        del(cpds_names_dict[cpd])

if(len(delete_cpds)>0):
    print("Removing "+str(len(delete_cpds))+" newly obsolete compounds")
    for cpd in delete_cpds:
        print("\t"+cpd)
    compounds_helper.saveCompounds(compounds_dict)
    compounds_helper.saveNames(cpds_names_dict)
    compounds_helper.saveAliases(cpds_aliases_dict)
