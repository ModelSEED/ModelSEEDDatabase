#!/usr/bin/env python
import os, sys, re, copy
from csv import DictReader
from collections import OrderedDict
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()

Names_Dict = CompoundsHelper.loadNames()
SearchNames_Dict = dict()
All_Names = dict()
New_Name_Count=dict()
for msid in sorted(Names_Dict.keys()):
    for name in Names_Dict[msid]:
        All_Names[name]=1

        searchname = CompoundsHelper.searchname(name)
        #Avoid redundancy where possible
        if(searchname not in SearchNames_Dict):
            SearchNames_Dict[searchname]=msid

Original_Alias_Dict=CompoundsHelper.loadMSAliases()
Source_Alias_Dict = dict()
All_Aliases=dict()
New_Alias_Count=dict()
for msid in Original_Alias_Dict.keys():
    for source in Original_Alias_Dict[msid].keys():
        if(source not in Source_Alias_Dict):
            Source_Alias_Dict[source]=dict()

        for alias in Original_Alias_Dict[msid][source]:
            if(alias not in All_Aliases):
                All_Aliases[alias]=dict()
            All_Aliases[alias][msid]=1

            if(alias not in Source_Alias_Dict[source]):
                Source_Alias_Dict[source][alias]=list()
            Source_Alias_Dict[source][alias].append(msid)

for alias in All_Aliases.keys():
    All_Aliases[alias]=sorted(All_Aliases[alias].keys())

Structures_Dict = CompoundsHelper.loadStructures(["InChI","SMILE"],["KEGG","MetaCyc"])
All_InChIs=dict()
All_Aliases_InChIs=dict()
for alias in Structures_Dict['InChI'].keys():
    if('Charged' in Structures_Dict['InChI'][alias]):
        for struct in Structures_Dict['InChI'][alias]['Charged'].keys():
            if(struct not in All_InChIs):
                All_InChIs[struct]=list()
            All_InChIs[struct].append(alias)

            if(alias not in All_Aliases_InChIs):
                All_Aliases_InChIs[alias]=list()
            All_Aliases_InChIs[alias].append(struct)
    elif('Original' in Structures_Dict['InChI'][alias]):
        for struct in Structures_Dict['InChI'][alias]['Original'].keys():
            if(struct not in All_InChIs):
                All_InChIs[struct]=list()
            All_InChIs[struct].append(alias)

            if(alias not in All_Aliases_InChIs):
                All_Aliases_InChIs[alias]=list()
            All_Aliases_InChIs[alias].append(struct)

All_SMILEs=dict()
All_Aliases_SMILEs=dict()
for alias in Structures_Dict['SMILE'].keys():
    if('Charged' in Structures_Dict['SMILE'][alias]):
        for struct in Structures_Dict['SMILE'][alias]['Charged'].keys():
            if(struct not in All_SMILEs):
                All_SMILEs[struct]=list()
            All_SMILEs[struct].append(alias)

            if(alias not in All_Aliases_SMILEs):
                All_Aliases_SMILEs[alias]=list()
            All_Aliases_SMILEs[alias].append(struct)
    elif('Original' in Structures_Dict['SMILE'][alias]):
        for struct in Structures_Dict['SMILE'][alias]['Original'].keys():
            if(struct not in All_SMILEs):
                All_SMILEs[struct]=list()
            All_SMILEs[struct].append(alias)

            if(alias not in All_Aliases_SMILEs):
                All_Aliases_SMILEs[alias]=list()
            All_Aliases_SMILEs[alias].append(struct)

#Find last identifier and increment
last_identifier = list(sorted(Compounds_Dict.keys()))[-1]
identifier_count = int(re.sub('^cpd','',last_identifier))

Biochem="KEGG"
Biochem_Root="../../Biochemistry/Aliases/Provenance/Primary_Databases/";

Default_Cpd = OrderedDict({ "id": "cpd00000","name": "null","abbreviation": "null","aliases": "null",
                             "formula": "null","mass": "10000000","charge": "0",
                             "deltag": "10000000","deltagerr": "10000000","pka": "","pkb": "",
                             "inchikey": "","smiles": "",
                             "is_cofactor": 0,"is_core": 0,"is_obsolete": 0,
                             "abstract_compound": "null","comprised_of": "null","linked_compound": "null",
                             "source": "" })

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
        if(cpd['ID'] in Source_Alias_Dict[Biochem]):
            matched_cpd = sorted(Source_Alias_Dict[Biochem][cpd['ID']])[0]
            matched_src="ID"

        #Then check that the Structure doesn't already exist, first as InChI, then as SMILE
        if(matched_cpd is None and cpd['InChI'] and cpd['InChI'] in All_InChIs):

            msids = dict()
            for alias in All_InChIs[cpd['InChI']]:

                #The structures are taken from their sources and the corresponding alias may not yet be registered
                if(alias not in All_Aliases):
                    continue

                for msid in All_Aliases[alias]:
                    msids[msid]=1

            msids=list(sorted(msids.keys()))
            if(len(msids)>0):
                matched_cpd=msids[0]
                matched_src='InChI'

        if(matched_cpd is None and cpd['SMILE'] and cpd['SMILE'] in All_SMILEs):

            msids = dict()
            for alias in All_SMILEs[cpd['SMILE']]:
                #The structures are taken from their sources and the corresponding alias may not yet be registered
                if(alias not in All_Aliases):
                    continue

                for msid in All_Aliases[alias]:
                    msids[msid]=1

            msids=list(sorted(msids.keys()))
            if(len(msids)>0):
                matched_cpd=msids[0]
                matched_src='SMILE'

        #Then check that the Name doesn't already exist
        if(matched_cpd is None):
            msids=dict()
            for name in cpd['NAMES'].split('|'):
                searchname = CompoundsHelper.searchname(name)
                if(searchname in SearchNames_Dict):
                    msids[SearchNames_Dict[searchname]]=1
            msids=list(sorted(msids.keys()))
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
                if(name not in All_Names):
                    #Possible for there to be no names in biochemistry?
                    if(matched_cpd not in Names_Dict):
                        Names_Dict[matched_cpd]=list()
                    Names_Dict[matched_cpd].append(name)
                    All_Names[name]=1
                    New_Name_Count[matched_cpd]=1

            #print warning if multiple structures
            if(cpd['InChI'] in All_InChIs):
                if(cpd['ID'] not in All_Aliases_InChIs or cpd['ID'] not in All_InChIs[cpd['InChI']]):
                    print("Warning: InChI structure for "+cpd['ID']+" assigned to different compounds: "+",".join(All_InChIs[cpd['InChI']]))

            #print warning if multiple structures
            if(cpd['SMILE'] in All_SMILEs):
                if(cpd['ID'] not in All_Aliases_SMILEs or cpd['ID'] not in All_SMILEs[cpd['SMILE']]):
                    print("Warning: SMILE structure for "+cpd['ID']+" assigned to different compounds: "+",".join(All_SMILEs[cpd['SMILE']]))
                
            #if matching structure or name, add ID to aliases
            if(matched_src != 'ID'):
                if(matched_cpd in Original_Alias_Dict and Biochem not in Original_Alias_Dict[matched_cpd]):
                    Original_Alias_Dict[matched_cpd][Biochem]=list()
                Original_Alias_Dict[matched_cpd][Biochem].append(cpd['ID'])
                New_Alias_Count[matched_cpd]=1

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

            #Add new identifier with KEGG ID as alias
            Original_Alias_Dict[new_cpd['id']]={Biochem:[cpd['ID']]}
            New_Alias_Count[new_cpd['id']]=1

            #Add new names
            #Names are saved separately as part of the aliases at the end of the script
            for name in cpd['NAMES'].split('|'):
                if(new_cpd['name']=='null'):
                    new_cpd['name']=name
                    new_cpd['abbreviation']=name

                if(name not in All_Names):
                    #Possible for there to be no names in biochemistry?
                    if(new_cpd['id'] not in Names_Dict):
                        Names_Dict[new_cpd['id']]=list()
                    Names_Dict[new_cpd['id']].append(name)
                    All_Names[name]=1
                    New_Name_Count[new_cpd['id']]=1

            #If no names at all
            if(new_cpd['name']=='null'):
                new_cpd['name']=cpd['ID']
                new_cpd['abbreviation']=cpd['ID']

            Compounds_Dict[new_cpd['id']]=new_cpd
            New_Cpd_Count[new_cpd['id']]=1

#Here, for matches, re-write names and aliases
print("Compounds matched via:")
for src in sorted(Matched_Cpd_Count.keys()):
    print("\t"+src+": "+str(len(Matched_Cpd_Count[src])))
print("Saving additional names for "+str(len(New_Name_Count))+" compounds")
CompoundsHelper.saveNames(Names_Dict)
print("Saving additional "+Biochem+" aliases for "+str(len(New_Alias_Count))+" compounds")
CompoundsHelper.saveAliases(Original_Alias_Dict)
print("Saving "+str(len(New_Cpd_Count.keys()))+" new compounds from "+Biochem)
CompoundsHelper.saveCompounds(Compounds_Dict)

#Scripts to run afterwards
#./Merge_Formulas.py
#./Update_Compound_Aliases.py
#../Structures/List_ModelSEED_Structures.py
#../Structures/Update_Compound_Structures.py
#./Update_Formula_Charge.py
