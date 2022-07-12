#!/usr/bin/env python
import os, sys, re, copy
import argparse, requests
from csv import DictReader
from collections import OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument('compounds_file', help="Compounds File")
parser.add_argument('github_user', help="GitHub username")
parser.add_argument("-s", dest='save_file', action='store_true')
parser.add_argument("-r", dest='report_file', action='store_true')
args = parser.parse_args()

if(os.path.isfile(args.compounds_file) is False):
    print("Cannot find file: "+args.compounds_file)
    sys.exit()

curation_source = args.github_user

url = f"https://api.github.com/users/{curation_source}"
r = requests.get(url.format(curation_source)).json()
if("login" in r and r["login"] == curation_source):
    print("Github user "+curation_source+" found")
if("message" in r and r["message"] == "Not Found"):
    print ("Github user "+curation_source+" does not exist.")
    sys.exit()

sys.path.append('../../../Libs/Python')
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

Structures_Dict = compounds_helper.loadStructures(["InChI","InChIKey","SMILE"],["KEGG","MetaCyc"])
all_structures=dict()
all_aliases_structures=dict()
for structure_type in ["InChI","InChIKey","SMILE"]:
    for alias in Structures_Dict[structure_type]:
        if('Charged' in Structures_Dict[structure_type][alias]):
            for structure in Structures_Dict[structure_type][alias]['Charged']:
                struct_list = [structure]
                if(structure_type == 'InChIKey'):
                    struct_list.append('-'.join(structure.split('-')[0:2]))
                    struct_list.append('-'.join(structure.split('-')[0:1]))

                for struct in struct_list:
                    if(struct not in all_structures):
                        all_structures[struct]=list()
                    all_structures[struct].append(alias)

                    if(alias not in all_aliases_structures):
                        all_aliases_structures[alias]=list()
                    all_aliases_structures[alias].append(struct)
                    
        elif('Original' in Structures_Dict[structure_type][alias]):
            for struct in Structures_Dict[structure_type][alias]['Original']:
                if(struct not in all_structures):
                    all_structures[struct]=list()
                all_structures[struct].append(alias)

                if(alias not in all_aliases_structures):
                    all_aliases_structures[alias]=list()
                all_aliases_structures[alias].append(struct)

#Find last identifier and increment
last_identifier = list(sorted(compounds_dict))[-1]
identifier_count = int(re.sub('^cpd','',last_identifier))

Default_Cpd = OrderedDict({ "id":"cpd00000","name":"null","abbreviation":"null","aliases":"null",
                             "formula":"null","mass":10000000,"charge":0,
                             "deltag":10000000.0,"deltagerr":10000000.0,"pka":"","pkb":"",
                             "inchikey":"","smiles":"",
                             "is_cofactor":0,"is_core":0,"is_obsolete":0,
                             "abstract_compound":"null","comprised_of":"null","linked_compound":"null",
                             "notes":[],"source":"" })
New_Cpd_Count=dict()
matched_cpds_dict=dict()
Headers=list()
with open(args.compounds_file) as fh:
    for line in fh.readlines():
        line=line.strip('\r\n')
        if(len(Headers)==0):
            Headers=line.split('\t')
            continue

        cpd=dict()
        array=line.split('\t',len(Headers))
        for i in range(len(Headers)):
            cpd[Headers[i].lower()]=array[i]

        (matched_cpd,matched_str,matched_src)=(None,None,None)

        #Check that the Structure doesn't already exist, first as InChI, then as InChiKey, then as SMILES
        if(matched_cpd is None):
            for structure_type in ['inchi','inchikey','smile','smiles']:
                if(structure_type in cpd):
                    structure_list = [cpd[structure_type]]
                    if(structure_type == 'inchikey'):
                        structure_list.append('-'.join(cpd[structure_type].split('-')[0:2]))
                        structure_list.append('-'.join(cpd[structure_type].split('-')[0:1]))

                    for structure in structure_list:
                        
                        if(matched_cpd is not None):
                            break

                        if(structure in all_structures):
                    
                            msids = list()
                            for alias in all_structures[structure]:

                                #The structures are taken from their sources and the corresponding alias may not yet be registered
                                if(alias not in all_aliases):
                                    continue

                                for msid in all_aliases[alias]:
                                    if(msid not in msids):
                                        msids.append(msid)

                            msids=list(sorted(msids))
                            if(len(msids)>0):
                                matched_cpd=msids[0]
                                matched_src=structure_type
                                matched_str=structure

        #Then check that the Name doesn't already exist
        elif(matched_cpd is None):
            msids=dict()
            for name in cpd['names'].split('|'):
                searchname = compounds_helper.searchname(name)
                if(searchname in searchnames_dict):
                    if(searchnames_dict[searchname] not in msids):
                        msids[searchnames_dict[searchname]]=name
            msids_list=list(sorted(msids))
            if(len(msids_list)>0):
                matched_cpd=msids_list[0]
                matched_src='NAMES'
                matched_str=msids[msids_list[0]]

        if(cpd['id'] not in matched_cpds_dict):
            matched_cpds_dict[cpd['id']]=list()
        match_dict = {'src':str(matched_src),'str':str(matched_str),'match':str(matched_cpd)}
        matched_cpds_dict[cpd['id']].append(match_dict)

        if(matched_cpd is not None):

            #Regardless of match-type, add new names
            #NB at this point, names shouldn't match _anything_ already in the database
            #Names are saved separately as part of the aliases at the end of the script
            for name in cpd['names'].split('|'):
                if(name not in all_names_dict):
                    #Possible for there to be no names in biochemistry?
                    if(matched_cpd not in names_dict):
                        names_dict[matched_cpd]=list()
                    names_dict[matched_cpd].append(name)
                    all_names_dict[name]=1
                    new_name_count[matched_cpd]=1

            #print warning if multiple structures
            if('inchi' in cpd and cpd['inchi'] in all_structures):
                if(cpd['id'] not in all_aliases_structures or cpd['id'] not in all_structures[cpd['inchi']]):
                    print("Warning: InChI structure for "+cpd['id']+" assigned to different compounds: "+",".join(all_structures[cpd['inchi']]))

            #print warning if multiple structures
            if('smile' in cpd and cpd['smile'] in all_structures):
                if(cpd['id'] not in all_aliases_structures or cpd['id'] not in all_structures[cpd['smile']]):
                    print("Warning: SMILE structure for "+cpd['id']+" assigned to different compounds: "+",".join(all_structures[cpd['smile']]))
                
            #if matching structure or name, add ID to aliases
            if(matched_src != 'ID'):
                if(matched_cpd not in original_alias_dict):
                    original_alias_dict[matched_cpd]=dict()
                if(matched_cpd in original_alias_dict and curation_source not in original_alias_dict[matched_cpd]):
                    original_alias_dict[matched_cpd][curation_source]=list()
                original_alias_dict[matched_cpd][curation_source].append(cpd['id'])
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
            if('mass' in cpd):
                new_cpd['mass']=float(cpd['mass'])
            if('charge' in cpd):
                new_cpd['charge']=int(cpd['charge'])
            if('formula' in cpd):
                new_cpd['formula']=cpd['formula']

            #Add new identifier with original ID as alias
            original_alias_dict[new_cpd['id']]={curation_source:[cpd['id']]}
            new_alias_count[new_cpd['id']]=1

            #Add new names
            #Names are saved separately as part of the aliases at the end of the script
            for name in cpd['names'].split('|'):
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
                new_cpd['name']=cpd['id']
                new_cpd['abbreviation']=cpd['id']

            #Add source type
            new_cpd['source']='Primary Database'
            compounds_dict[new_cpd['id']]=new_cpd
            New_Cpd_Count[new_cpd['id']]=1

#Here, for matches, re-write names and aliases
print("Compounds matched via:")
matched_src_dict=dict()
for oc in sorted(matched_cpds_dict):
    for match_dict in matched_cpds_dict[oc]:
        if(match_dict['src'] not in matched_src_dict):
            matched_src_dict[match_dict['src']]=list()
        matched_src_dict[match_dict['src']].append(match_dict['match'])

for src in matched_src_dict:
    if(src != 'None'):
        match_dict = dict()
        for match in matched_src_dict[src]:
            match_dict[match]=1
        print("\t"+src+" "+str(len(matched_src_dict[src]))+" matched to "+str(len(match_dict.keys())))
    else:
        print("\t"+src+" "+str(len(matched_src_dict[src])))
        
if(args.report_file is True):
    file_stub = '.'.join(args.compounds_file.split('.')[0:-1])
    report_file = file_stub+'.rpt'
    print("Saving report to file: "+report_file)
    with open(report_file,'w') as rfh:
        for oc in sorted(matched_cpds_dict):
            for match_dict in matched_cpds_dict[oc]:
                rfh.write(oc+'\t'+match_dict['match']+'\t'+match_dict['src']+'\t'+match_dict['str']+'\n')
        
if(args.save_file is True):
    print("Saving additional names for "+str(len(new_name_count))+" compounds")
    compounds_helper.saveNames(names_dict)
    print("Saving additional "+curation_source+" aliases for "+str(len(new_alias_count))+" compounds")
    compounds_helper.saveAliases(original_alias_dict)
    print("Saving "+str(len(New_Cpd_Count))+" new compounds from "+curation_source)
    compounds_helper.saveCompounds(compounds_dict)

#Scripts to run afterwards
#./Merge_Formulas.py
#./Update_Compound_Aliases.py
#../Structures/List_ModelSEED_Structures.py
#../Structures/Update_Compound_Structures_Formulas_Charge.py
#./Rebalance_Reactions.py
