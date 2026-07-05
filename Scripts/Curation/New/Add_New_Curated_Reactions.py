#!/usr/bin/env python
import os, sys, re, copy
import argparse, requests
from collections import OrderedDict

def list_of_strings(arg):
    return arg.split(',')

parser = argparse.ArgumentParser()
parser.add_argument('reactions_file', help="Reactions File")
parser.add_argument('cpd_database', help="Biochemistry database of origin for compounds", type=list_of_strings)
parser.add_argument('github_user', help="GitHub username")
parser.add_argument("-r", dest='report_file', action='store_true')
parser.add_argument("-s", dest='save_file', action='store_true')
args = parser.parse_args()

if(os.path.isfile(args.reactions_file) is False):
    print("Cannot find file: "+args.reactions_file)
    sys.exit()
    
curation_source = args.github_user

url = f"https://api.github.com/users/{curation_source}"
r = requests.get(url.format(curation_source)).json()
if("login" in r and r["login"].lower() == curation_source.lower()):
    print("Github user "+curation_source+" found")
else:
    print ("Github user "+curation_source+" does not exist.")
    sys.exit()

sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

compounds_helper = Compounds()
compounds_alias_dict=compounds_helper.loadMSAliases()
source_alias_dict = dict()
for msid in compounds_alias_dict:
    for source in compounds_alias_dict[msid]:
        if(source not in source_alias_dict):
            source_alias_dict[source]=dict()
        for alias in compounds_alias_dict[msid][source]:
            if(alias not in source_alias_dict[source]):
                source_alias_dict[source][alias]=list()
            source_alias_dict[source][alias].append(msid)

#Check compound source
for db in args.cpd_database:
    if(db != "ModelSEED" and db not in source_alias_dict):
        print("Alias ("+db+") for source of compounds is not recognized")
        sys.exit()

compounds_dict = compounds_helper.loadCompounds()
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
reactions_codes = reactions_helper.generateCodes(reactions_dict)

Default_Rxn = OrderedDict({"id":"cpd00001","abbreviation":"null","name":"null","aliases":"null",
                           "code":"null","stoichiometry":"null","equation":"null","definition":"null",
                           "reversibility":"=","direction":"=","deltag":"10000000","deltagerr":"10000000",
                           "status":"NB","is_obsolete":0,"is_transport":0,
                           "abstract_reaction":"null","pathways":"null","ec_numbers":"null",
                           "compound_ids":"null","linked_reaction":"null","notes":[],"source":""})

original_rxn_alias_dict=reactions_helper.loadMSAliases()
new_alias_count=dict()

original_name_dict = reactions_helper.loadNames()
All_Names = dict()
New_Name_Count=dict()
for msid in sorted(original_name_dict):
    for name in original_name_dict[msid]:
        All_Names[name]=1

original_ec_dict = reactions_helper.loadECs()
All_ECs = dict()
New_EC_Count=dict()
for msid in sorted(original_ec_dict):
    for name in original_ec_dict[msid]:
        All_ECs[name]=1

#Find last identifier and increment
last_identifier = list(sorted(reactions_dict))[-1]
identifier_count = int(re.sub('^rxn','',last_identifier))

#If a reaction, after removing cpd redundancy, is empty
#We use a placeholder
Empty_Rxn_ID="rxn14003"

New_Rxn_Count=dict()
Headers=list()
missing_cpds=dict()
matched_rxns_dict=dict()
with open(args.reactions_file) as fh:
    for line in fh.readlines():
        line=line.strip()
        if(len(Headers)==0):
            Headers=line.split('\t')
            continue

        rxn=dict()
        array=line.split('\t',len(Headers))
        for i in range(len(Headers)):
            rxn[Headers[i].lower()]=array[i]

        #Retrieve identifiers from within equation
        #Split based on whitespace, and remove compartment index
        original_cpd_array=rxn['equation'].split(' ')

        new_cpd_array=list()
        for i in range(len(original_cpd_array)):
            if(re.search(r'^\([\d\.]+\)$',original_cpd_array[i])):
                continue

            if(original_cpd_array[i] == '+' ):
                continue

            if(re.search(r'^<?[-=]+>?$',original_cpd_array[i])):
                continue

            if(re.search(r'\[[01]\]$',original_cpd_array[i])):
                new_cpd_array.append(re.sub(r'\[[01]\]$','',original_cpd_array[i]))
            else:
                new_cpd_array.append(original_cpd_array[i])

        all_matched=True
        eqn_missing_cpds=list()
        for new_cpd in new_cpd_array:
            
            #if source is ModelSEED
            msid = ""
            if("ModelSEED" in args.cpd_database and new_cpd in compounds_dict):
                msid = new_cpd
            else:
                for db in args.cpd_database:
                    if(db in source_alias_dict and new_cpd in source_alias_dict[db]):
                        msid = sorted(source_alias_dict[db][new_cpd])[0]

            if(msid != ""):

                #set boundary
                bound_msid=msid+"["
                bound_cpd=new_cpd+"["
                esc_cpd = re.escape(bound_cpd)
                
                eqn_array = rxn['equation'].split(" ")
                new_eqn_array = list()
                for entry in eqn_array:
                    entry = re.sub(esc_cpd,bound_msid,entry)
                    new_eqn_array.append(entry)
                rxn['equation'] = " ".join(new_eqn_array)
            else:
                missing_cpds[new_cpd]=1
                eqn_missing_cpds.append(new_cpd)
                all_matched=False

        if(all_matched is False):
            print("Warning: missing "+",".join(args.cpd_database)+" identifiers for reaction "+rxn['id']+": "+"|".join(eqn_missing_cpds))
            continue
        
        rxn_cpds_array = reactions_helper.parseEquation(rxn['equation'])
        adjusted=False
        new_rxn_cpds_array = reactions_helper.removeCpdRedundancy(rxn_cpds_array)
        if(len(new_rxn_cpds_array)!=len(rxn_cpds_array)):
            adjusted=True
        
        rxn_code = reactions_helper.generateCode(new_rxn_cpds_array)
        matched_rxn=None        
        if(len(new_rxn_cpds_array)==0):
            matched_rxn = Empty_Rxn_ID
        else:
            if(rxn_code in reactions_codes):
                matched_rxn = sorted(list(reactions_codes[rxn_code]))[0]

        #Because we adjust for water a posterior
        #We need to include water when matching codes, in case
        if(matched_rxn is None):

            #Find statuses that only have water imbalance
            new_status = reactions_helper.balanceReaction(new_rxn_cpds_array)
            if(new_status == "MI:H:2/O:1" or new_status == "MI:H:-2/O:-1"):
                Water_Adjustment = 1
                if("-1" in new_status):
                    Water_Adjustment = -1

                #Adjust for water
                reactions_helper.adjustCompound(new_rxn_cpds_array,"cpd00001",float(Water_Adjustment))
                rxn_code = reactions_helper.generateCode(new_rxn_cpds_array)
                if(rxn_code in reactions_codes):
                    matched_rxn = sorted(list(reactions_codes[rxn_code]))[0]

        if(rxn['id'] not in matched_rxns_dict):
            matched_rxns_dict[rxn['id']]=list()
        matched_rxns_dict[rxn['id']].append(matched_rxn)

        if(matched_rxn is not None):
            #Add Names, EC and Alias
            #Regardless of match-type, add new names
            #NB at this point, names shouldn't match _anything_ already in the database
            #Names are saved separately as part of the aliases at the end of the script
            if('names' in rxn):
                for name in rxn['names'].split('|'):
                    if(name not in All_Names):
                        #Possible for there to be no names in biochemistry?
                        if(matched_rxn not in original_name_dict):
                            original_name_dict[matched_rxn]=list()
                        original_name_dict[matched_rxn].append(name)
                        All_Names[name]=1
                        New_Name_Count[matched_rxn]=1

            if('ecs' in rxn):
                for ec in rxn['ecs'].split('|'):
                    if(ec not in All_ECs):
                        #Possible for there to be no ecs in biochemistry?
                        if(matched_rxn not in original_ec_dict):
                            original_ec_dict[matched_rxn]=list()
                        original_ec_dict[matched_rxn].append(ec)
                        All_ECs[ec]=1
                        New_EC_Count[matched_rxn]=1

            #Add ID to aliases if the match is with a different reaction
            if(matched_rxn not in original_rxn_alias_dict):
                original_rxn_alias_dict[matched_rxn]=dict()
            if(matched_rxn in original_rxn_alias_dict and curation_source not in original_rxn_alias_dict[matched_rxn]):
                original_rxn_alias_dict[matched_rxn][curation_source]=list()
            if(rxn['id'] not in original_rxn_alias_dict[matched_rxn][curation_source]):
                original_rxn_alias_dict[matched_rxn][curation_source].append(rxn['id'])
                new_alias_count[matched_rxn]=1

            #Update source type
            reactions_dict[matched_rxn]['source']='User'

        elif(args.save_file is True):
            
            #New Reaction!
            #Generate new identifier
            identifier_count+=1
            new_identifier = 'rxn'+str(identifier_count)

            new_rxn = copy.deepcopy(Default_Rxn)
            new_rxn['id']=new_identifier

            #Add new identifier with compound source as alias
            original_rxn_alias_dict[new_rxn['id']]={curation_source:[rxn['id']]}
            new_alias_count[new_rxn['id']]=1

            #Add new names
            #Names are saved separately as part of the aliases at the end of the script
            if('names' in rxn):
                for name in rxn['names'].split('|'):
                    if(new_rxn['name']=='null'):
                        new_rxn['name']=name
                        new_rxn['abbreviation']=name

                    if(name not in All_Names):
                        #Possible for there to be no names in biochemistry?
                        if(new_rxn['id'] not in original_name_dict):
                            original_name_dict[new_rxn['id']]=list()
                        original_name_dict[new_rxn['id']].append(name)
                        All_Names[name]=1
                        New_Name_Count[new_rxn['id']]=1

            #If no names at all
            if(new_rxn['name']=='null'):
                new_rxn['name']=rxn['id']
                new_rxn['abbreviation']=rxn['id']

            #ECs
            if('ecs' in rxn):
                for ec in rxn['ecs'].split('|'):
                    if(ec not in All_ECs):
                        #Possible for there to be no ecs in biochemistry?
                        if(new_rxn['id'] not in original_ec_dict):
                            original_ec_dict[new_rxn['id']]=list()
                        original_ec_dict[new_rxn['id']].append(ec)
                        All_ECs[ec]=1
                        New_EC_Count[new_rxn['id']]=1

            #Add source type
            new_rxn['source']='User'

            reactions_dict[new_rxn['id']]=new_rxn
            New_Rxn_Count[new_rxn['id']]=1

            #Rebuild key fields for reaction using parsed equation
            reactions_helper.rebuildReaction(reactions_dict[new_rxn['id']],new_rxn_cpds_array)

            #Finally, because several new reactions may share equations
            if(rxn_code not in reactions_codes):
                reactions_codes[rxn_code]=dict()
            reactions_codes[rxn_code][new_rxn['id']]=1

if(len(missing_cpds)>0):
    print("Warning, these compounds in reactions were missing in the database")
    print("\t"+"|".join(sorted(missing_cpds))+"\n")

matched_dict={'match':[],None:[]}
for oc in sorted(matched_rxns_dict):
    for matched_rxn in matched_rxns_dict[oc]:
        if(matched_rxn is not None):
            matched_dict['match'].append(matched_rxn)
        else:
            matched_dict[None].append(matched_rxn)

msid_match_dict=dict()
for mc in matched_dict['match']:
    msid_match_dict[mc]=1
print("Reactions matched: "+str(len(matched_dict['match']))+" matched to "+str(len(msid_match_dict.keys())))
print("\t"+str(len(matched_dict[None]))+" not matched to any ModelSEED reactions")

if(args.report_file is True):
    file_stub = '.'.join(args.reactions_file.split('.')[0:-1])
    report_file = file_stub+'.rpt'
    print("Saving report to file: "+report_file)
    with open(report_file,'w') as rfh:
        for oc in sorted(matched_rxns_dict):
            for match in matched_rxns_dict[oc]:
                rfh.write(oc+'\t'+str(match)+'\n')
                
if(args.save_file is True):
    #Here, for matches, re-write names, ecs and aliases
    print("Saving additional ECs for "+str(len(New_EC_Count))+" reactions")
    reactions_helper.saveECs(original_ec_dict)
    print("Saving additional names for "+str(len(New_Name_Count))+" reactions")
    reactions_helper.saveNames(original_name_dict)
    print("Saving additional "+curation_source+" aliases for "+str(len(new_alias_count))+" reactions")
    reactions_helper.saveAliases(original_rxn_alias_dict)
    print("Saving "+str(len(New_Rxn_Count))+" new reactions from "+curation_source)
    reactions_helper.saveReactions(reactions_dict)

#Scripts to run afterwards
#../../Biochemistry/Refresh/Rebalance_Reactions.py (very important)
#../../Biochemistry/Refresh/Adjust_Reaction_Protons.py
#../../Biochemistry/Refresh/Adjust_Reaction_Water.py
#../../Biochemistry/Refresh/Merge_Reactions.py (merges may happen because of water)
