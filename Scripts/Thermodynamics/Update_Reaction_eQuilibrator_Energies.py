#!/usr/bin/env python
import os,sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions

label="eQuilibrator"
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

        eq_reactions[array[0]]=[float("{0:.2f}".format(float(array[1]))),float("{0:.2f}".format(float(array[2])))]

file_handle.close()

# print(len(eq_reactions))
# 13, 874 ModelSEED Reactions

for rxn in sorted (reactions_dict.keys()):

    if(rxn not in eq_reactions):
        #NB There are a number of reactions for which there are structures available
        #for every eQuilibrator record (and labeled EQC in reaction notes)
        #but, for whatever reason, we couldn't retrieve an estimated energy for
        #the reaction. This is likely because we couldn't retrieve an estimated
        #energy for every compound with a structure in eQuilibrator
        continue

    #Here we establish an arbitrary threshold of 100 for the error, if the error
    #is too big, we shouldn't use it, but for storing the database we keep it now
    if(float(eq_reactions[rxn][1]) > 100):
        pass

    # values always saved as list of energy and error
    if(not isinstance(reactions_dict[rxn]['thermodynamics'],dict)):
        reactions_dict[rxn]['thermodynamics'] = dict()
    if(label not in reactions_dict[rxn]['thermodynamics']):
        reactions_dict[rxn]['thermodynamics'][label]=list()
    reactions_dict[rxn]['thermodynamics'][label]=eq_reactions[rxn]

print("Saving reactions")
reactions_helper.saveReactions(reactions_dict)