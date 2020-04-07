#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

# We actually don't want obsolete reactions and compounds in our database
# So we're striving to remove any 'new' ones that are obsolete
# Any information attached to them should be associated with their linked counterpart

last_rxn_str='rxn40000'
last_rxn_int=int(last_rxn_str[3:])

delete_rxns=list()
for rxn in reactions_dict:
    rxn_int = int(rxn[3:])
    if(rxn_int > last_rxn_int and reactions_dict[rxn]['is_obsolete']):
        delete_rxns.append(rxn)

for rxn in delete_rxns:

    # I need to update the linked_reaction field to make sure that the removed reaction is deleted
    for lnkd_rxn in reactions_dict[rxn]['linked_reaction']:
        reactions_dict[lnkd_rxn]['linked_reaction'].remove(rxn)

    del reactions_dict[rxn]

if(len(delete_rxns)>0):
    print("Removing "+str(len(delete_rxns))+" newly obsolete reactions")
    reactions_helper.saveReactions(reactions_dict)
