#!/usr/bin/env python
import os, sys, re, copy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('reactions_file', help="Reactions File")
parser.add_argument('compound_source', help="ModelSEED alias or GitHub username")
parser.add_argument("-s", dest='save_file', action='store_true')
args = parser.parse_args()



if(os.path.isfile(args.reactions_file) is False):
    print("Cannot find file: "+reactions_file)
    sys.exit()

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

compounds_helper = Compounds()
Compounds_Alias_Dict=compounds_helper.loadMSAliases()
Source_Alias_Dict = dict()
for msid in Compounds_Alias_Dict:
    for source in Compounds_Alias_Dict[msid]:
        if(source not in Source_Alias_Dict):
            Source_Alias_Dict[source]=dict()
        for alias in Compounds_Alias_Dict[msid][source]:
            if(alias not in Source_Alias_Dict[source]):
                Source_Alias_Dict[source][alias]=list()
            Source_Alias_Dict[source][alias].append(msid)

#Check compound source
if(args.compound_source != "ModelSEED" and args.compound_source not in Source_Alias_Dict):
    print("Alias for source of compounds is not recognized")
    sys.exit()

compounds_dict = compounds_helper.loadCompounds()
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
reactions_codes = reactions_helper.generateCodes(reactions_dict)

Default_Rxn = {"id":"cpd00001","name":"null","abbreviation":"null","aliases":[],
               "code":"null","stoichiometry":"null","equation":"null","definition":"null",
               "reversibility":"=","direction":"=","deltag":"10000000","deltagerr":"10000000",
               "status":"NB","is_obsolete":0,"is_transport":0,
               "abstract_reaction":"null","pathways":"null","ec_numbers":"null",
               "compound_ids":"null","linked_reaction":"null","notes":[],"source":""}

original_rxn_alias_dict=reactions_helper.loadMSAliases()
New_Alias_Count=dict()

Names_Dict = reactions_helper.loadNames()
All_Names = dict()
New_Name_Count=dict()
for msid in sorted(Names_Dict):
    for name in Names_Dict[msid]:
        All_Names[name]=1

ECs_Dict = reactions_helper.loadECs()
All_ECs = dict()
New_EC_Count=dict()
for msid in sorted(ECs_Dict):
    for name in ECs_Dict[msid]:
        All_ECs[name]=1

#Find last identifier and increment
last_identifier = list(sorted(reactions_dict))[-1]
identifier_count = int(re.sub('^rxn','',last_identifier))

#If a reaction, after removing cpd redundancy, is empty
#We use a placeholder
Empty_Rxn_ID="rxn14003"

New_Rxn_Count=dict()
Headers=list()
rxns=list()
missing_cpds=dict()
with open(args.reactions_file) as fh:
    for line in fh.readlines():
        line=line.strip()
        if(len(Headers)==0):
            Headers=line.split('\t')
            continue

        rxn=dict()
        array=line.split('\t',len(Headers))
        for i in range(len(Headers)):
            rxn[Headers[i].upper()]=array[i]

        #Retrieve identifiers from within equation
        #Split based on whitespace, and remove compartment index
        original_cpd_array=rxn['EQUATION'].split(' ')

        new_cpd_array=list()
        for i in range(len(original_cpd_array)):
            if(re.search('^\([\d\.]+\)$',original_cpd_array[i])):
                continue

            if(original_cpd_array[i] == '+' ):
                continue

            if(re.search('^<?[-=]+>?$',original_cpd_array[i])):
                continue

            if(re.search('\[[01]\]$',original_cpd_array[i])):
                new_cpd_array.append(re.sub('\[[01]\]$','',original_cpd_array[i]))
            else:
                new_cpd_array.append(original_cpd_array[i])

        all_matched=True
        eqn_missing_cpds=list()
        for new_cpd in new_cpd_array:
            
            #if source is ModelSEED
            msid = ""
            if(args.compound_source == 'ModelSEED' and new_cpd in compounds_dict):
                msid = new_cpd
            elif(args.compound_source in Source_Alias_Dict and new_cpd in Source_Alias_Dict[args.compound_source]):
                msid = sorted(Source_Alias_Dict[args.compound_source][new_cpd])[0]

            if(msid != ""):

                #set boundary
                bound_msid=msid+"["
                bound_cpd=new_cpd+"["
                esc_cpd = re.escape(bound_cpd)
                
                eqn_array = rxn['EQUATION'].split(" ")
                new_eqn_array = list()
                for entry in eqn_array:
                    entry = re.sub(esc_cpd,bound_msid,entry)
                    new_eqn_array.append(entry)
                rxn['EQUATION'] = " ".join(new_eqn_array)
            else:
                missing_cpds[new_cpd]=1
                eqn_missing_cpds.append(new_cpd)
                all_matched=False

        if(all_matched is False):
            print("Warning: missing "+args.compound_source+" identifiers for reaction "+rxn['ID']+": "+"|".join(eqn_missing_cpds))
            continue
        
        rxn_cpds_array = reactions_helper.parseEquation(rxn['EQUATION'])
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

        if(matched_rxn is not None):
            print("Matched reaction:\t"+rxn['ID']+"\t"+matched_rxn)
            #Add Names, EC and Alias
            #Regardless of match-type, add new names
            #NB at this point, names shouldn't match _anything_ already in the database
            #Names are saved separately as part of the aliases at the end of the script
            if('NAMES' in rxn):
                for name in rxn['NAMES'].split('|'):
                    if(name not in All_Names):
                        #Possible for there to be no names in biochemistry?
                        if(matched_rxn not in Names_Dict):
                            Names_Dict[matched_rxn]=list()
                        Names_Dict[matched_rxn].append(name)
                        All_Names[name]=1
                        New_Name_Count[matched_rxn]=1

            if('ECS' in rxn):
                for ec in rxn['ECS'].split('|'):
                    if(ec not in All_ECs):
                        #Possible for there to be no ecs in biochemistry?
                        if(matched_rxn not in ECs_Dict):
                            ECs_Dict[matched_rxn]=list()
                        ECs_Dict[matched_rxn].append(ec)
                        All_ECs[ec]=1
                        New_EC_Count[matched_rxn]=1

            #Add ID to aliases if the match is with a different reaction
            if(matched_rxn not in original_rxn_alias_dict):
                original_rxn_alias_dict[matched_rxn]=dict()
            if(matched_rxn in original_rxn_alias_dict and args.compound_source not in original_rxn_alias_dict[matched_rxn]):
                original_rxn_alias_dict[matched_rxn][args.compound_source]=list()
            if(rxn['ID'] not in original_rxn_alias_dict[matched_rxn][args.compound_source]):
                original_rxn_alias_dict[matched_rxn][args.compound_source].append(rxn['ID'])
                New_Alias_Count[matched_rxn]=1

            #Update source type
            reactions_dict[matched_rxn]['source']=args.compound_source

        else:
            
            #New Reaction!
            #Generate new identifier
            identifier_count+=1
            new_identifier = 'rxn'+str(identifier_count)

            new_rxn = copy.deepcopy(Default_Rxn)
            new_rxn['id']=new_identifier

            #Add new identifier with compound source as alias
            original_rxn_alias_dict[new_rxn['id']]={args.compound_source:[rxn['ID']]}
            New_Alias_Count[new_rxn['id']]=1

            #Add new names
            #Names are saved separately as part of the aliases at the end of the script
            if('NAMES' in rxn):
                for name in rxn['NAMES'].split('|'):
                    if(new_rxn['name']=='null'):
                        new_rxn['name']=name
                        new_rxn['abbreviation']=name

                    if(name not in All_Names):
                        #Possible for there to be no names in biochemistry?
                        if(new_rxn['id'] not in Names_Dict):
                            Names_Dict[new_rxn['id']]=list()
                        Names_Dict[new_rxn['id']].append(name)
                        All_Names[name]=1
                        New_Name_Count[new_rxn['id']]=1

            #If no names at all
            if(new_rxn['name']=='null'):
                new_rxn['name']=rxn['ID']
                new_rxn['abbreviation']=rxn['ID']

            #ECs
            if('ECS' in rxn):
                for ec in rxn['ECS'].split('|'):
                    if(ec not in All_ECs):
                        #Possible for there to be no ecs in biochemistry?
                        if(new_rxn['id'] not in ECs_Dict):
                            ECs_Dict[new_rxn['id']]=list()
                        ECs_Dict[new_rxn['id']].append(ec)
                        All_ECs[ec]=1
                        New_EC_Count[new_rxn['id']]=1

            #Add source type
            new_rxn['source']=args.compound_source

            reactions_dict[new_rxn['id']]=new_rxn
            New_Rxn_Count[new_rxn['id']]=1

            #Rebuild key fields for reaction using parsed equation
            stoichiometry=reactions_helper.buildStoich(new_rxn_cpds_array)
            reactions_helper.rebuildReaction(reactions_dict[new_rxn['id']],stoichiometry)

            #Finally, because several new reactions may share equations
            if(rxn_code not in reactions_codes):
                reactions_codes[rxn_code]=dict()
            reactions_codes[rxn_code][new_rxn['id']]=1
#            print(new_rxn['id'],rxn_code,rxn['ID'])

if(len(missing_cpds)>0):
    print("Missing Compounds: "+"|".join(sorted(missing_cpds)))

#Here, for matches, re-write names, ecs and aliases
if(len(New_EC_Count)>0):
    print("Saving additional ECs for "+str(len(New_EC_Count))+" reactions")
    if(args.save_file is True):
        reactions_helper.saveECs(ECs_Dict)

if(len(New_Name_Count)>0):
    print("Saving additional names for "+str(len(New_Name_Count))+" reactions")
    if(args.save_file is True):
        reactions_helper.saveNames(Names_Dict)

if(len(New_Alias_Count)>0):
    print("Saving additional "+args.compound_source+" aliases for "+str(len(New_Alias_Count))+" reactions")
    if(args.save_file is True):
        reactions_helper.saveAliases(original_rxn_alias_dict)

if(len(New_Rxn_Count)>0):
    print("Saving "+str(len(New_Rxn_Count))+" new reactions from "+args.compound_source)
    if(args.save_file is True):
        reactions_helper.saveReactions(reactions_dict)

#Scripts to run afterwards
#../../Biochemistry/Refresh/Rebalance_Reactions.py (very important)
#../../Biochemistry/Refresh/Adjust_Reaction_Protons.py
#../../Biochemistry/Refresh/Adjust_Reaction_Water.py
#../../Biochemistry/Refresh/Merge_Reactions.py (merges may happen because of water)
