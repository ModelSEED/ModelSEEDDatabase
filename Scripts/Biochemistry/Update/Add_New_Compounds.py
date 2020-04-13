#!/usr/bin/env python
import os, sys, re, copy
from csv import DictReader
from collections import OrderedDict
temp=list();
header=True;

Biochem="MetaCyc"

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

names_dict = compounds_helper.loadNames()
searchnames_dict = dict()
all_names_dict = dict()
new_name_count = dict()
for msid in sorted(names_dict):
    for name in names_dict[msid]:
        all_names_dict[name]=1

        searchname = compounds_helper.searchname(name)
        #Avoid redundancy where possible
        if(searchname not in searchnames_dict):
            searchnames_dict[searchname]=msid

original_alias_dict=compounds_helper.loadMSAliases()
source_alias_dict = dict()
all_aliases = dict()
new_alias_count = dict()
for msid in original_alias_dict:
    for source in original_alias_dict[msid]:
        if(source not in source_alias_dict):
            source_alias_dict[source]=dict()

        for alias in original_alias_dict[msid][source]:
            if(alias not in all_aliases):
                all_aliases[alias]=dict()
            all_aliases[alias][msid]=1

            if(alias not in source_alias_dict[source]):
                source_alias_dict[source][alias]=list()
            source_alias_dict[source][alias].append(msid)

for alias in all_aliases:
    all_aliases[alias]=sorted(all_aliases[alias])

Structures_Dict = compounds_helper.loadStructures(["InChI","SMILE"],["KEGG","MetaCyc"])
all_inchis=dict()
all_aliases_InChIs=dict()
for alias in Structures_Dict['InChI']:
    if('Charged' in Structures_Dict['InChI'][alias]):
        for struct in Structures_Dict['InChI'][alias]['Charged']:
            if(struct not in all_inchis):
                all_inchis[struct]=list()
            all_inchis[struct].append(alias)

            if(alias not in all_aliases_InChIs):
                all_aliases_InChIs[alias]=list()
            all_aliases_InChIs[alias].append(struct)
    elif('Original' in Structures_Dict['InChI'][alias]):
        for struct in Structures_Dict['InChI'][alias]['Original']:
            if(struct not in all_inchis):
                all_inchis[struct]=list()
            all_inchis[struct].append(alias)

            if(alias not in all_aliases_InChIs):
                all_aliases_InChIs[alias]=list()
            all_aliases_InChIs[alias].append(struct)

all_smiles=dict()
all_aliases_SMILEs=dict()
for alias in Structures_Dict['SMILE']:
    if('Charged' in Structures_Dict['SMILE'][alias]):
        for struct in Structures_Dict['SMILE'][alias]['Charged']:
            if(struct not in all_smiles):
                all_smiles[struct]=list()
            all_smiles[struct].append(alias)

            if(alias not in all_aliases_SMILEs):
                all_aliases_SMILEs[alias]=list()
            all_aliases_SMILEs[alias].append(struct)
    elif('Original' in Structures_Dict['SMILE'][alias]):
        for struct in Structures_Dict['SMILE'][alias]['Original']:
            if(struct not in all_smiles):
                all_smiles[struct]=list()
            all_smiles[struct].append(alias)

            if(alias not in all_aliases_SMILEs):
                all_aliases_SMILEs[alias]=list()
            all_aliases_SMILEs[alias].append(struct)

#Find last identifier and increment
last_identifier = list(sorted(compounds_dict))[-1]
identifier_count = int(re.sub('^cpd','',last_identifier))

Biochem_Root="../../Biochemistry/Aliases/Provenance/Primary_Databases/";

Default_Cpd = OrderedDict({ "id":"cpd00000","name":"null","abbreviation":"null","aliases":"null",
                             "formula":"null","mass":"10000000","charge":"0",
                             "deltag":"10000000","deltagerr":"10000000","pka":"","pkb":"",
                             "inchikey":"","smiles":"",
                             "is_cofactor":0,"is_core":0,"is_obsolete":0,
                             "abstract_compound":"null","comprised_of":"null","linked_compound":"null",
                             "source":"" })

Matched_Cpd_Count=dict()
New_Cpd_Count=dict()
Headers=list()
cpds=list()
with open(Biochem_Root+Biochem+"_Compounds.tbl") as fh:
    for line in fh.readlines():
        line=line.strip()
        if(len(Headers)==0):
            Headers=line.split('\t')
            continue

        cpd=dict()
        array=line.split('\t',len(Headers))
        for i in range(len(Headers)):
            cpd[Headers[i]]=array[i]

        (matched_cpd,matched_src)=(None,None)

        #First check that the Alias doesn't already exist
        if(cpd['ID'] in source_alias_dict[Biochem]):
            matched_cpd = sorted(source_alias_dict[Biochem][cpd['ID']])[0]
            matched_src="ID"

        #Then check that the Structure doesn't already exist, first as InChI, then as SMILE
        if(matched_cpd is None and cpd['InChI'] and cpd['InChI'] in all_inchis):

            msids = dict()
            for alias in all_inchis[cpd['InChI']]:

                #The structures are taken from their sources and the corresponding alias may not yet be registered
                if(alias not in all_aliases):
                    continue

                for msid in all_aliases[alias]:
                    msids[msid]=1

            msids=list(sorted(msids))
            if(len(msids)>0):
                matched_cpd=msids[0]
                matched_src='InChI'

        elif(matched_cpd is None and cpd['SMILE'] and cpd['SMILE'] in all_smiles):

            msids = dict()
            for alias in all_smiles[cpd['SMILE']]:
                #The structures are taken from their sources and the corresponding alias may not yet be registered
                if(alias not in all_aliases):
                    continue

                for msid in all_aliases[alias]:
                    msids[msid]=1

            msids=list(sorted(msids))
            if(len(msids)>0):
                matched_cpd=msids[0]
                matched_src='SMILE'

        #Then check that the Name doesn't already exist
        elif(matched_cpd is None):
            msids=dict()
            for name in cpd['NAMES'].split('|'):
                searchname = compounds_helper.searchname(name)
                if(searchname in searchnames_dict):
                    msids[searchnames_dict[searchname]]=1
            msids=list(sorted(msids))
            if(len(msids)>0):
                matched_cpd=msids[0]
                matched_src='NAMES'

        if(matched_cpd is not None):
            
            if(matched_src not in Matched_Cpd_Count):
                Matched_Cpd_Count[matched_src]=list()
            Matched_Cpd_Count[matched_src].append(matched_cpd)

            #Regardless of match-type, add new names
            #NB at this point, names shouldn't match _anything_ already in the database
            #Names are saved separately as part of the aliases at the end of the script
            for name in cpd['NAMES'].split('|'):
                if(name not in all_names_dict):
                    #Possible for there to be no names in biochemistry?
                    if(matched_cpd not in names_dict):
                        names_dict[matched_cpd]=list()
                    names_dict[matched_cpd].append(name)
                    all_names_dict[name]=1
                    new_name_count[matched_cpd]=1

            #print warning if multiple structures
            if(cpd['InChI'] in all_inchis):
                if(cpd['ID'] not in all_aliases_InChIs or cpd['ID'] not in all_inchis[cpd['InChI']]):
                    print("Warning: InChI structure for "+cpd['ID']+" assigned to different compounds: "+",".join(all_inchis[cpd['InChI']]))

            #print warning if multiple structures
            if(cpd['SMILE'] in all_smiles):
                if(cpd['ID'] not in all_aliases_SMILEs or cpd['ID'] not in all_smiles[cpd['SMILE']]):
                    print("Warning: SMILE structure for "+cpd['ID']+" assigned to different compounds: "+",".join(all_smiles[cpd['SMILE']]))
                
            #if matching structure or name, add ID to aliases
            if(matched_src != 'ID'):
                if(matched_cpd not in original_alias_dict):
                    original_alias_dict[matched_cpd]={Biochem:list()}
                if(matched_cpd in original_alias_dict and Biochem not in original_alias_dict[matched_cpd]):
                    original_alias_dict[matched_cpd][Biochem]=list()
                original_alias_dict[matched_cpd][Biochem].append(cpd['ID'])
                new_alias_count[matched_cpd]=1

            #Update source type
            compounds_dict[matched_cpd]['source']='Primary Database'

        else:

            #New Compound!
            #Generate new identifier
            identifier_count+=1
            new_identifier = 'cpd'+str(identifier_count)

            new_cpd = copy.deepcopy(Default_Cpd)
            new_cpd['id']=new_identifier
            new_cpd['mass']=cpd['MASS']
            new_cpd['charge']=cpd['CHARGE']
            new_cpd['formula']=cpd['FORMULA']

            #Add new identifier with original ID as alias
            original_alias_dict[new_cpd['id']]={Biochem:[cpd['ID']]}
            new_alias_count[new_cpd['id']]=1

            #Add new names
            #Names are saved separately as part of the aliases at the end of the script
            for name in cpd['NAMES'].split('|'):
                if(new_cpd['name']=='null'):
                    new_cpd['name']=name
                    new_cpd['abbreviation']=name

                if(name not in all_names_dict):
                    #Possible for there to be no names in biochemistry?
                    if(new_cpd['id'] not in names_dict):
                        names_dict[new_cpd['id']]=list()
                    names_dict[new_cpd['id']].append(name)
                    all_names_dict[name]=1
                    new_name_count[new_cpd['id']]=1

            #If no names at all
            if(new_cpd['name']=='null'):
                new_cpd['name']=cpd['ID']
                new_cpd['abbreviation']=cpd['ID']

            #Add source type
            new_cpd['source']='Primary Database'
            compounds_dict[new_cpd['id']]=new_cpd
            New_Cpd_Count[new_cpd['id']]=1

#Here, for matches, re-write names and aliases
print("Compounds matched via:")
for src in sorted(Matched_Cpd_Count):
    print("\t"+src+": "+str(len(Matched_Cpd_Count[src])))
print("Saving additional names for "+str(len(new_name_count))+" compounds")
compounds_helper.saveNames(names_dict)
print("Saving additional "+Biochem+" aliases for "+str(len(new_alias_count))+" compounds")
compounds_helper.saveAliases(original_alias_dict)
print("Saving "+str(len(New_Cpd_Count))+" new compounds from "+Biochem)
compounds_helper.saveCompounds(compounds_dict)

#Scripts to run afterwards
#./Merge_Formulas.py
#./Update_Compound_Aliases.py
#../Structures/List_ModelSEED_Structures.py
#../Structures/Update_Compound_Structures_Formulas_Charge.py
#./Rebalance_Reactions.py
