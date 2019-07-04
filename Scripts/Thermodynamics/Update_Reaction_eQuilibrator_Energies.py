#!/usr/bin/env python
import os,sys
from BiochemPy import Reactions

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

############################################################################
##
## We apply/overwrite with eQuilibrator Energies
##
############################################################################

thermodynamics_root=os.path.dirname(__file__)+"/../../Biochemistry/Thermodynamics/"
file_name=thermodynamics_root+'eQuilibrator/MetaNetX_Reaction_Energies.tbl'
eq_reactions=dict()
with open(file_name) as file_handle:
    for line in file_handle.readlines():
        line = line.strip()
        array= line.split('\t')

        if('energy' in array[1] or array[1] == 'nan'):
            continue

        eq_reactions[array[0]]={'dg':"{0:.2f}".format(float(array[1])),'dge':"{0:.2f}".format(float(array[2]))}

file_handle.close()

# print(len(eq_reactions))
# 13, 874 ModelSEED Reactions

file_handle = open('Reactions_GroupFormation_eQuilibrator_Comparison.txt', 'w')
file_handle.write('ID\tGF\tEQ\n')
for rxn in sorted (reactions_dict.keys()):

    rxn_gf='nan'
    rxn_eq='nan'

    if("GFC" in reactions_dict[rxn]['notes']):
        rxn_gf='|'.join([str(reactions_dict[rxn]['deltag']),str(reactions_dict[rxn]['deltagerr'])])

    if(rxn in eq_reactions):
        rxn_eq='|'.join([eq_reactions[rxn]['dg'],eq_reactions[rxn]['dge']])

    #Write the values to file. I'm not making exception, but I'm
    #including this `if` statement as a comment to remind me of ways in which
    #they can't be directly compared
    #if(reactions_dict[rxn]['deltag']!=10000000 and rxn_gf != 'nan' and rxn_eq != 'nan'):
    file_handle.write('\t'.join([rxn,rxn_gf,rxn_eq])+'\n')

    #Having printed to file, we skip where there is no eQuilibrator estimate
    if(rxn not in eq_reactions):
        #NB There are a number of reactions for which there are structures available
        #for every eQuilibrator record (and labeled EQC in reaction notes)
        #but, for whatever reason, we couldn't retrieve an estimated energy for
        #the reaction. This is likely because we couldn't retrieve an estimated
        #energy for every compound with a structure in eQuilibrator
        continue

    #Here we establish an arbitrary threshold of 100 for the error, if the error
    #is too big, we don't use it

    if(float(eq_reactions[rxn]['dge']) > 100):
        continue

    notes_list=reactions_dict[rxn]['notes']
    if(not isinstance(notes_list,list)):
        notes_list=list()

    reactions_dict[rxn]['deltag']=float(eq_reactions[rxn]['dg'])
    reactions_dict[rxn]['deltagerr']=float(eq_reactions[rxn]['dge'])
    if('EQU' not in notes_list):
        notes_list.append('EQU')
    reactions_dict[rxn]['notes']=notes_list

file_handle.close()
print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)
