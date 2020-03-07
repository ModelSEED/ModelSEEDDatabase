#!/usr/bin/env python
import os, sys
from BiochemPy import Compounds, Reactions, InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

parents_children=dict()
children_parents=dict()
inchi_exception="WQZGKKKJIJFFOK"
cpd_exception="cpd00027"
for cpd in compounds_dict:
    cpd_obj = compounds_dict[cpd]
    if( not isinstance(cpd_obj['ontology'],dict) ):
        continue

    if( inchi_exception not in cpd_obj['inchikey'] ):
        continue

    print(cpd_obj['id'],cpd_obj['inchikey'])
    if('parent_class' in cpd_obj['ontology']):
        for cpd in cpd_obj['ontology']['parent_class'].split(";"):
            if(cpd_obj['id'] not in children_parents):
                children_parents[cpd_obj['id']]=dict()
            children_parents[cpd_obj['id']][cpd]=1
            if(cpd not in parents_children):
                parents_children[cpd]=dict()
            parents_children[cpd][cpd_obj['id']]=1

#    print( cpd_obj['id'], cpd_obj['ontology'], isinstance(cpd_obj['ontology'],dict) )
print(parents_children)

#in first pass, find all reactions that contain either children or parents
reactions_to_use={'parent':[],'child':[],'mixed':[]}
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

for rxn in reactions_dict:
    rxn_obj = reactions_dict[rxn]
    for cpd in rxn_obj['compound_ids'].split(';'):

        if( cpd != cpd_exception ):
            continue

        if(cpd in parents_children):
            reactions_to_use['parent'].append(rxn)
        if(cpd in children_parents):
            reactions_to_use['child'].append(rxn)

#First take all child reactions, and generate a set where all relationships are parents
for rxn in reactions_to_use['child']:
    rxn_obj = reactions_dict[rxn]
#    print(rxn_obj['id'])
    if(rxn in reactions_to_use['child']):
        #indicates that a reaction has a mixture of parent and children
        pass
    
    rxn_cpds_array = reactions_helper.parseStoich(rxn_obj['stoichiometry'])
    for rgt in rxn_cpds_array:
        if(rgt['compound'] in children_parents):
            print(rgt['compound'],children_parents[rgt['compound']])

    
    #new_rxn_cpds_array = reactions_helper.removeCpdRedundancy(rxn_cpds_array)
    break

for rxn in reactions_to_use['parent']:
    rxn_obj = reactions_dict[rxn]
    print(rxn_obj['id'],rxn_obj['definition'])

sys.exit()
compounds_helper = Compounds()
structures_dict = compounds_helper.loadStructures(["InChI"],["ModelSEED"])
compounds_dict = compounds_helper.loadCompounds()

#Notes for Sets of Stereoisomers
#38 WQZGKKKJIJFFOK hexoses 6C (1 generic missing: L-Mannose)
#20 SRBFZHDQGSBBOR aldopentoses (1 generic missing: D-Lyxose)
#20 SHZGCJCMOBCMKK deoxyhexoses
#14 OVRNDRQMDRJTHS aminohexoses
#14 HXXFSFRBOHSIMQ hexose-1-P
#14 AEMOLEFTQBMNLQ hexuronic acid
#13 GUBGYTABKSRVRQ hexose dimer
#11 RGHNJXZEOKUKBD hexonic acid
#11 RFSUNEUAIZKAJO hex-ketose
#11 DLRVVLDZNNYCBX hexose dimer
#10 NBSCHQHZLSJFNQ hexose-6-P
#10 KRCZYMFUWVJCLI monoterpene
#10 HMFHBZSHGGEWLO pentose
#10 GZCGUPFRVQAUEE more hexoses (open form)
#10 FBPFZTCFMRRESA alcohol hexose

substructures = {"WQZGKKKJIJFFOK":6,"SRBFZHDQGSBBOR":5,"SHZGCJCMOBCMKK":6}

for substructure in substructures:
    n_carbons = substructures[substructure]
    print("==============================")
    groups = dict()
    stereos = dict()
    for cpd in compounds_dict:
        cpd_obj = compounds_dict[cpd]
        if(cpd_obj['is_obsolete']==0 and substructure in cpd_obj['inchikey']):

            #Parse out InChI formula and layers
            (inchi_formula,inchi_layers) = InChIs.parse(list(structures_dict[cpd]['InChI'])[0])

            #We're using stereochemistry so skip if we don't have the layer
            if(inchi_layers['t'] == ""):
                continue

            all_stereos = inchi_layers['t'].split(',')
            stereos[cpd_obj['id']]=all_stereos
        
#           print(cpd_obj['id'],cpd_obj['name'],inchi_layers)

            short_stereos = list()
            for number in range(n_carbons-3,n_carbons):
                for stereo in all_stereos:
                    if(str(number) in stereo):
                        short_stereos.append(stereo)
            short_stereos.append(inchi_layers['m'])
            if("".join(short_stereos) not in groups):
                groups["".join(short_stereos)]=list()
            groups["".join(short_stereos)].append(cpd_obj)

    for group in groups:

        #Only deal with sets of grouped compounds greater than 1
        if(len(groups[group])<2):
            continue

        parents = list()
        children = list()
        for obj in sorted(groups[group], key = lambda x: x['id']):

            #Ignore if not comparing stereoisomers
            if(len(stereos[obj['id']])==0):
                continue

            #There's a single incorrect instance of a D/L Glucose having a '?' in the first position
            if('?' in stereos[obj['id']][0]):
                continue

            #Sugars have '?' at last position if alpha/beta generic
            if('?' in stereos[obj['id']][-1]):
                parents.append(obj['id'])
            else:
                children.append(obj['id'])

            print("\t".join([group,obj['id'], obj['name'], str(stereos[obj['id']][0]), str(stereos[obj['id']][-1])]))

        print("\t".join([substructure,group,str(parents),str(children),"\n"]))
        if(len(parents)>0 and len(children)>0):

            #Add children class for parent
            for parent in parents:
                if(compounds_dict[parent]['ontology'] == "null"):
                    compounds_dict[parent]['ontology']=dict()

                compounds_dict[parent]['ontology']['child_class']=";".join(children)

            #Add parent class for children
            for child in children:
                if(compounds_dict[child]['ontology'] == "null"):
                    compounds_dict[child]['ontology']=dict()

                compounds_dict[child]['ontology']['parent_class']=";".join(parents)

        elif(group != "" and len(groups[group])>0):
            print("No parents and or children found for group '"+group+"'")

print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
