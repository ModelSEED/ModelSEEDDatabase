#!/usr/bin/env python
import sys
sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions, Compounds
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

cpd_file="../../../Biochemistry/Aliases/Provenance/Rhea/ChEBI_ID_Name_InChIKey.rpt"
name_matches=list()
with open(cpd_file) as fh:
    for line in fh.readlines():
        line=line.strip('\r\n')
        array=line.split('\t')
        if(array[2] == "name"):
            name_matches.append(array[1])

old_file="../../../Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt"
old_match=dict()
with open(old_file) as fh:
    for line in fh.readlines():
        line=line.strip('\r\n')
        array=line.split('\t')
        if(array[2] == "rhea"):
            if(array[1] not in old_match):
                old_match[array[1]]=list()
            old_match[array[1]].append(array[0])

new_file="../../../Biochemistry/Aliases/Provenance/Rhea/Rhea_reactions.rpt"
new_match=dict()
with open(new_file) as fh:
    for line in fh.readlines():
        line=line.strip('\r\n')
        array=line.split('\t')
        if(array[1] == 'None'):
            continue
        if(array[0] not in new_match):
            new_match[array[0]]=list()
        new_match[array[0]].append(array[1])

global_cpds=dict()
for new_rhea in new_match:
    #perfect match
    if(new_rhea in old_match):
        if(new_match[new_rhea] == old_match[new_rhea]):
            pass
        else:
#            print(new_rhea,new_match[new_rhea],old_match[new_rhea])

            new_cpds=list()
            for new_ms in new_match[new_rhea]:
                for new_cpd in reactions_dict[new_ms]['compound_ids'].split(';'):
                    if(new_cpd not in new_cpds):
                        new_cpds.append(new_cpd)
                        
            old_cpds=list()
            for old_ms in old_match[new_rhea]:
                for old_cpd in reactions_dict[old_ms]['compound_ids'].split(';'):
                    if(old_cpd not in old_cpds):
                        old_cpds.append(old_cpd)

            for new_cpd in new_cpds:
                if(new_cpd not in old_cpds):
                    if(new_cpd not in global_cpds):
                        global_cpds[new_cpd]=list()
                    global_cpds[new_cpd].append(new_rhea)
                    
            for old_cpd in old_cpds:
                if(old_cpd not in new_cpds):
                    if(old_cpd not in global_cpds):
                        global_cpds[old_cpd]=list()
                    global_cpds[old_cpd].append(new_rhea)
            pass
    else:
        pass

#for rhea in global_cpds['cpd28218']:
#    if(rhea in global_cpds['cpd00004']):
#    if(rhea not in global_cpds['cpd00087']):
#        print(rhea,new_match[rhea],old_match[rhea])

skip = list()
# Spot-checked some NAD-OR-NADP reactions, these are changed for actual NAD or NADP reactions
skip.append('cpd27638')
skip.append('cpd27640')

# Spot-checked some NADP/NADPH reactions, these were correct
skip.append('cpd00005')
skip.append('cpd00006')

# Spot-checked some NAD/NADH reactions, these were correct
skip.append('cpd00003')
skip.append('cpd00004')

# Not checking protons and water and oxygen
skip.append('cpd00067')
skip.append('cpd00001')
skip.append('cpd00007')

# Spot-checked THF and THF-Glu-N generic
skip.append('cpd00087')
skip.append('cpd28218')

# Spot-checked D-Glucose and b-D-Glucose
skip.append('cpd00027')
skip.append('cpd00190')

# Spot-checked Acetate - need to remove "ester" from names and merge to cpd11966
skip.append('cpd00029')

# Checked Menaquinol
skip.append('cpd11451')

# Spot-checked D-Fructose
skip.append('cpd00082')

# Spot-checked Reduced flavin
skip.append('cpd21035')

# Spot-checked Oxo-pentoate and 2-Dehydro-3-deoxy-L-pentonate
skip.append('cpd32214')
skip.append('cpd00461')
skip.append('cpd00520')

# Spot-checked Pyruvate
skip.append('cpd00020')

# Spot-checked RCN
skip.append('cpd00542')

# Spot-checked 1-Pyrroline-2-carboxylate
skip.append('cpd02235')

# Spot-checked N-acetylneuraminate (several stereoisomers)
skip.append('cpd00232')

# Spot-checked N-acetyl-D-Mannosamine
skip.append('cpd00492')

# Spot-checked 5-Dehydro-D-fructose
skip.append('cpd00234')

# Spot-checked Copalyl diphosphate
skip.append('cpd03629')

# Spot-checked Galactose
skip.append('cpd00108')

# Spot-checked Phenylcarbinol
skip.append('cpd00435')

# Spot-checked 1,2-Diacyl-sn-glycerol
skip.append('cpd11423')

# Generics
# cpd11610 and cpd11609 and cpd22290! and cpd26978! are generic
# cpd00049 and cpd26675 are generic carboxylic acid
# cpd11611 and cpd22234 are generic acyl-CoA
# cpd00057 and cpd22318 are generic alcohol
skip.append('cpd11610')
skip.append('cpd11609')
skip.append('cpd22290')
skip.append('cpd26978')
skip.append('cpd00049')
skip.append('cpd26675')
skip.append('cpd11611')
skip.append('cpd22234')
skip.append('cpd00057')
skip.append('cpd22318')

for cpd in sorted(global_cpds.keys(), key=lambda vl: len(global_cpds[vl])):
    if(cpd not in skip and cpd in name_matches):
        if(len(global_cpds[cpd])>0):
            print("\n==============================================================================")
            print(cpd,"\t",compounds_dict[cpd]['name'],"\t",len(global_cpds[cpd]),"\t",global_cpds[cpd])
            for rhea in global_cpds[cpd]:
                if(rhea in old_match):
                    print("\tOld: ",rhea,"\t",",".join(old_match[rhea]),"\t",reactions_dict[old_match[rhea][0]]['definition'])
                if(rhea in new_match):
                    print("\tNew: ",rhea,"\t",",".join(new_match[rhea]),"\t",reactions_dict[new_match[rhea][0]]['definition'])

#        if(array[0] in old_match and array[1] in old_match[array[0]]):
#            print(array[0],array[1])
#            match_count+=1
#        if(array[0] in old_match and array[1] not in old_match[array[0]]):
#            print(array[0],array[1],old_match[array[0]])
# print(len(old_match.keys()),match_count)
# match_count=0
