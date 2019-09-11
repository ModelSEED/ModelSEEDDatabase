#!/usr/bin/env python
import os, sys, re, copy
from csv import DictReader
from collections import OrderedDict
temp=list();
header=True;

if(len(sys.argv)<5):
    print("Not enough arguments!")
    print("./Add_New_Compounds.py <file> <biochemistry> <primary identifiers:(0|1)> <source> <prefix>")
    sys.exit()

Biochem_File=sys.argv[1]
Biochem_Dir=Biochem_File.split('/')[0]
Primary_Biochem=sys.argv[2]
Primary_IDs=int(sys.argv[3])
Biochem_Source=sys.argv[4]
Biochem_Source="Published Model"
ID_Prefix="cpd"
if(len(sys.argv)==6):
    ID_Prefix=sys.argv[5]+ID_Prefix

Biochem=Primary_Biochem
#If it isn't the actual biochemistry database, then use proper name
if(Biochem not in Biochem_File):
    array=Biochem_Dir.split('_')
    Biochem=array[0]
    if(array[1].isdigit() is False):
        Biochem='_'.join([Biochem,array[1]])

#Find curated compound mappings
import os.path
curated_mappings=dict()
if(os.path.exists(Biochem_Dir+'/Curated_Compound_Mappings.tsv')):
    with open(Biochem_Dir+'/Curated_Compound_Mappings.tsv') as fh:
        for line in fh.readlines():
            line=line.strip()
            (original_cpd,ms_cpd)=line.split('\t')
            curated_mappings[original_cpd]=ms_cpd
    fh.close()

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
identifier_count=0
identifiers = list(sorted(filter(lambda x:ID_Prefix in x, compounds_dict)))
if(len(identifiers)>0):
    last_identifier = identifiers[-1]
    identifier_count = int(re.sub('^\w*cpd','',last_identifier))

Default_Cpd = OrderedDict({ "id":"cpd00000","name":"null","abbreviation":"null","aliases":"null",
                             "formula":"null","mass":"10000000","charge":"0",
                             "deltag":"10000000","deltagerr":"10000000","pka":"","pkb":"",
                             "inchikey":"","smiles":"",
                             "is_cofactor":0,"is_core":0,"is_obsolete":0,
                             "abstract_compound":"null","comprised_of":"null","linked_compound":"null",
                             "source":"", "ontology":"null", "notes":"null" })

Matched_Cpd_Count=dict()
New_Cpd_Count=dict()
Headers=list()
cpds=list()
cpd_integration_report=dict()
with open(Biochem_File) as fh:
    for line in fh.readlines():
        line=line.strip('\r\n')

        if(len(Headers)==0):
            Headers=line.split('\t')
            continue

        cpd=dict()
        array=line.split('\t',len(Headers))
        for i in range(len(Headers)):
            cpd[Headers[i]]=array[i]
            
        cpd_integration_report[cpd['ID']]=array

        (matched_cpd,matched_src)=(None,None)

        #First check that the compound hasn't been manually curated
        if(cpd['ID'] in curated_mappings):
            matched_cpd=curated_mappings[cpd['ID']]
            matched_src="MANUAL"

        id_to_check = cpd['ID']
        if(matched_cpd is None and Primary_Biochem != "InChI"):
            #Check to see if using primary id or ids in another column
            if(Primary_IDs==0 and Primary_Biochem in cpd and cpd[Primary_Biochem] != ''):
                id_to_check=cpd[Primary_Biochem]

            #First check that the Alias doesn't already exist
            if(id_to_check in source_alias_dict[Primary_Biochem]):
                matched_cpd = sorted(source_alias_dict[Primary_Biochem][id_to_check])[0]
                if(Primary_IDs==1):
                    matched_src="ID"
                else:
                    matched_src=Primary_Biochem

        #Then check that the Structure doesn't already exist, first as InChI, then as SMILE
        if(matched_cpd is None and 'InChI' in cpd and cpd['InChI'] and cpd['InChI'] in all_inchis):

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

        elif(matched_cpd is None and 'SMILE' in cpd and cpd['SMILE'] and cpd['SMILE'] in all_smiles):

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
            for name in re.split('[|;]',cpd['NAMES']):
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
            for name in re.split('[|;]',cpd['NAMES']):
                if(name not in all_names_dict):
                    #Possible for there to be no names in biochemistry?
                    if(matched_cpd not in names_dict):
                        names_dict[matched_cpd]=list()
                    names_dict[matched_cpd].append(name)
                    all_names_dict[name]=1
                    new_name_count[matched_cpd]=1

                    #here, check and add searchname too, if its new
                    searchname = compounds_helper.searchname(name)
                    if(searchname not in searchnames_dict):
                        searchnames_dict[searchname]=matched_cpd

            #print warning if multiple structures
            if('InChI' in cpd and cpd['InChI'] in all_inchis):
                if(id_to_check not in all_aliases_InChIs or id_to_check not in all_inchis[cpd['InChI']]):
                    print("Warning: InChI structure for "+id_to_check+" assigned to different compounds: "+",".join(all_inchis[cpd['InChI']]))

            #print warning if multiple structures
            if('SMILE' in cpd and cpd['SMILE'] in all_smiles):
                if(id_to_check not in all_aliases_SMILEs or id_to_check not in all_smiles[cpd['SMILE']]):
                    print("Warning: SMILE structure for "+id_to_check+" assigned to different compounds: "+",".join(all_smiles[cpd['SMILE']]))
                
            #if matching structure or name, add ID to aliases
            if(matched_cpd not in original_alias_dict):
                original_alias_dict[matched_cpd]={Biochem:list()}
            if(Biochem not in original_alias_dict[matched_cpd]):
                original_alias_dict[matched_cpd][Biochem]=list()
            if(cpd['ID'] not in original_alias_dict[matched_cpd][Biochem]):
                original_alias_dict[matched_cpd][Biochem].append(cpd['ID'])
                new_alias_count[matched_cpd]=1

            #Update source type
            if(Biochem_Source=='Primary Database'):
                compounds_dict[matched_cpd]['source']=Biochem_Source
            elif(Biochem_Source=='Secondary Database' and compounds_dict[matched_cpd]['source'] != 'Primary Database'):
                compounds_dict[matched_cpd]['source']=Biochem_Source
            elif(Biochem_Source=='Published Model' and 'Database' not in compounds_dict[matched_cpd]['source']):
                compounds_dict[matched_cpd]['source']=Biochem_Source

            #Matched
            if(matched_cpd in New_Cpd_Count):
                #Can't be a match if its with a new compound (duplicate matches are OK)
                cpd_integration_report[cpd['ID']].append('N')
            else:
                cpd_integration_report[cpd['ID']].append('Y')
            cpd_integration_report[cpd['ID']].append(matched_cpd)
            cpd_integration_report[cpd['ID']].append(compounds_dict[matched_cpd]['name'])

        else:

            #New Compound!
            #Generate new identifier
            identifier_count+=1
            new_identifier = ID_Prefix+str(identifier_count).zfill(5)

            new_cpd = copy.deepcopy(Default_Cpd)
            new_cpd['id']=new_identifier
            if('MASS' in cpd):
                new_cpd['mass']=cpd['MASS']
            if('CHARGE' in cpd):
                new_cpd['charge']=cpd['CHARGE']
            if('FORMULA' in cpd):
                new_cpd['formula']=cpd['FORMULA']

            #Add new identifier with original ID as alias
            original_alias_dict[new_cpd['id']]={Biochem:[cpd['ID']]}
            new_alias_count[new_cpd['id']]=1

            #Add new names
            #Names are saved separately as part of the aliases at the end of the script
            for name in re.split('[|;]',cpd['NAMES']):
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

                    #here, check and add searchname too, if its new
                    searchname = compounds_helper.searchname(name)
                    if(searchname not in searchnames_dict):
                        searchnames_dict[searchname]=new_cpd['id']

            #If no names at all
            if(new_cpd['name']=='null'):
                new_cpd['name']=cpd['ID']
                new_cpd['abbreviation']=cpd['ID']

            #Add source type
            new_cpd['source']=Biochem_Source
            compounds_dict[new_cpd['id']]=new_cpd
            New_Cpd_Count[new_cpd['id']]=1

            #Matched
            cpd_integration_report[cpd['ID']].append('N')
            cpd_integration_report[cpd['ID']].append(new_cpd['id'])
            cpd_integration_report[cpd['ID']].append('')

with open(Biochem_Dir+'/'+Biochem+'_Compound_Integration_Report.txt','w') as fh:
    fh.write('\t'.join(Headers)+'\t'+'\t'.join(['MATCH','MODELSEED','MS NAMES'])+'\n')
    for cpd in sorted(cpd_integration_report):
        fh.write('\t'.join(cpd_integration_report[cpd])+'\n')
fh.close()

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
