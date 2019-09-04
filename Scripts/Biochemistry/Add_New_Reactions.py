#!/usr/bin/env python
import os, sys, re, copy
from csv import DictReader
from collections import OrderedDict
temp=list();
header=True;

if(len(sys.argv)<3):
    print("Not enough arguments!")
    print("./Add_New_Reactions.py <file> <biochemistry> <source> <prefix>")
    sys.exit()

Biochem_File=sys.argv[1]
Biochem_Dir=Biochem_File.split('/')[0]
Biochem_Source=sys.argv[2]
Biochem_Source="Published Model"
ID_Prefix=""
if(len(sys.argv)==4):
    ID_Prefix=sys.argv[3]

Biochem=""
array=Biochem_Dir.split('_')
Biochem=array[0]
if(array[1].isdigit() is False):
    Biochem='_'.join([Biochem,array[1]])

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
reactions_codes = reactions_helper.generateCodes(reactions_dict)

Default_Rxn = {"id":"cpd00001","name":"null","abbreviation":"null","aliases":"null",
               "code":"null","stoichiometry":"null","equation":"null","definition":"null",
               "reversibility":"=","direction":"=","deltag":"10000000","deltagerr":"10000000",
               "status":"NB","is_obsolete":0,"is_transport":0,
               "abstract_reaction":"null","pathways":"null","ec_numbers":"null",
               "compound_ids":"null","linked_reaction":"null","notes":"null","source":""}

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

Original_Alias_Dict=reactions_helper.loadMSAliases()
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
mismatched_cpds=dict()
rxn_integration_report=dict()
with open(Biochem_File) as fh:
    for line in fh.readlines():
        line=line.strip('\r\n')
        if(len(Headers)==0):
            Headers=line.split('\t')
            continue

        rxn=dict()
        array=line.split('\t',len(Headers))
        for i in range(len(Headers)):
            rxn[Headers[i]]=array[i]

        rxn_integration_report[rxn['ID']]=array[0:2]

        is_transport='N'
        if('[0]' in rxn['EQUATION'] and '[1]' in rxn['EQUATION']):
            is_transport='Y'

        #To replace identifiers within equation
        #They are retrieved and then have the compartment index removed
        original_cpd_array=rxn['EQUATION'].split(' ')
        original_stripped_cpd_array=list()
        for i in range(len(original_cpd_array)):
            if(re.search('\[[01]\]$',original_cpd_array[i])):
                original_stripped_cpd_array.append(re.sub('\[[01]\]$','',original_cpd_array[i]))

        #For rebuilding equation array that's useful for integration
        rxn['REPORT_EQUATION']=rxn['EQUATION']

        original_matched_cpds=list()
        for original_cpd in original_stripped_cpd_array:
            if(original_cpd in Source_Alias_Dict[Biochem]):
                msid = sorted(Source_Alias_Dict[Biochem][original_cpd])[0]

                #set boundary for regular expression
                bound_msid=msid+"["
                bound_cpd=original_cpd+"["
                esc_cpd = re.escape(bound_cpd)
                esc_cpd="^"+esc_cpd

                eqn_array = rxn['EQUATION'].split(" ")
                match_eqn_array = list()
                for entry in eqn_array:
                    match_entry = re.sub(esc_cpd,bound_msid,entry)
                    match_eqn_array.append(match_entry)
                rxn['EQUATION'] = " ".join(match_eqn_array)

                #if original ms compound, update report
                if(msid[0:3] == 'cpd'):
                    original_matched_cpds.append(msid)

                    bound_msid=original_cpd+"("+msid+")["
                    eqn_array = rxn['REPORT_EQUATION'].split(" ")
                    report_eqn_array = list()
                    for entry in eqn_array:
                        report_entry = re.sub(esc_cpd,bound_msid,entry)
                        report_eqn_array.append(report_entry)
                    rxn['REPORT_EQUATION'] = " ".join(report_eqn_array)

        rxn_cpds_array = reactions_helper.parseEquation(rxn['EQUATION'])
        adjusted='N'
        new_rxn_cpds_array = reactions_helper.removeCpdRedundancy(rxn_cpds_array)
        if(len(new_rxn_cpds_array)!=len(rxn_cpds_array)):
            adjusted='Y'

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
                    adjusted='W'

        rxn_integration_report[rxn['ID']].append(rxn['REPORT_EQUATION']) #equation
        rxn_integration_report[rxn['ID']].append(is_transport) #transport
        rxn_integration_report[rxn['ID']].append(adjusted) #adjustments made before match

        if(matched_rxn is not None and matched_rxn not in New_Rxn_Count):
            rxn_integration_report[rxn['ID']].append('Y') #match
        else:
            rxn_integration_report[rxn['ID']].append('N') #match

            eqn_array = rxn['REPORT_EQUATION'].split(" ")
            for entry in eqn_array:
                if(re.search('\[[01]\]$',entry) and 'cpd' not in entry):
                    if(entry not in mismatched_cpds):
                        mismatched_cpds[entry]=list()
                    mismatched_cpds[entry].append(rxn['ID'])

        if(len(original_matched_cpds)==len(original_stripped_cpd_array)):
            rxn_integration_report[rxn['ID']].append('Y') #complete
        else:
            rxn_integration_report[rxn['ID']].append('N') #complete

        if(matched_rxn is not None):
            #Add Names, EC and Alias
            #Regardless of match-type, add new names
            #NB at this point, names shouldn't match _anything_ already in the database
            #Names are saved separately as part of the aliases at the end of the script
            for name in rxn['NAMES'].split('|'):
                if(name not in All_Names):
                    #Possible for there to be no names in biochemistry?
                    if(matched_rxn not in Names_Dict):
                        Names_Dict[matched_rxn]=list()
                    Names_Dict[matched_rxn].append(name)
                    All_Names[name]=1
                    New_Name_Count[matched_rxn]=1

            if('ECs' in rxn):
                for ec in rxn['ECs'].split('|'):
                    if(ec not in All_ECs):
                        #Possible for there to be no ecs in biochemistry?
                        if(matched_rxn not in ECs_Dict):
                            ECs_Dict[matched_rxn]=list()
                        ECs_Dict[matched_rxn].append(ec)
                        All_ECs[ec]=1
                        New_EC_Count[matched_rxn]=1

            #Add ID to aliases if the match is with a different reaction
            if(matched_rxn not in Original_Alias_Dict):
                Original_Alias_Dict[matched_rxn]=dict()
            if(matched_rxn in Original_Alias_Dict and Biochem not in Original_Alias_Dict[matched_rxn]):
                Original_Alias_Dict[matched_rxn][Biochem]=list()
            if(rxn['ID'] not in Original_Alias_Dict[matched_rxn][Biochem]):
                Original_Alias_Dict[matched_rxn][Biochem].append(rxn['ID'])
                New_Alias_Count[matched_rxn]=1

            #Update source type
            if(Biochem_Source=='Primary Database'):
                reactions_dict[matched_rxn]['source']=Biochem_Source
            elif(Biochem_Source=='Secondary Database' and reactions_dict[matched_rxn]['source'] != 'Primary Database'):
                reactions_dict[matched_rxn]['source']=Biochem_Source
            elif(Biochem_Source=='Published Model' and 'Database' not in reactions_dict[matched_rxn]['source']):
                reactions_dict[matched_rxn]['source']=Biochem_Source

            #Matched
            rxn_integration_report[rxn['ID']].append(matched_rxn)
            rxn_integration_report[rxn['ID']].append(reactions_dict[matched_rxn]['definition'])

        else:
            
            #New Reaction!
            #Generate new identifier
            identifier_count+=1
            new_identifier = ID_Prefix+'rxn'+str(identifier_count)

            new_rxn = copy.deepcopy(Default_Rxn)
            new_rxn['id']=new_identifier

            #Add new identifier with KEGG ID as alias
            Original_Alias_Dict[new_rxn['id']]={Biochem:[rxn['ID']]}
            New_Alias_Count[new_rxn['id']]=1

            #Add new names
            #Names are saved separately as part of the aliases at the end of the script
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
            if('ECs' in rxn):
                for ec in rxn['ECs'].split('|'):
                    if(ec not in All_ECs):
                        #Possible for there to be no ecs in biochemistry?
                        if(new_rxn['id'] not in ECs_Dict):
                            ECs_Dict[new_rxn['id']]=list()
                        ECs_Dict[new_rxn['id']].append(ec)
                        All_ECs[ec]=1
                        New_EC_Count[new_rxn['id']]=1

            #Add source type
            new_rxn['source']=Biochem_Source

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

            #Matched
            rxn_integration_report[rxn['ID']].append(new_rxn['id'])
            rxn_integration_report[rxn['ID']].append('')

with open(Biochem_Dir+'/'+Biochem+'_Reaction_Integration_Report.txt','w') as fh:
    fh.write('\t'.join(Headers)+'\t'+'\t'.join(['TRANSPORT','ADJUSTED','MATCH','COMPLETE','MODELSEED','MS DEFINITION'])+'\n')
    for cpd in sorted(rxn_integration_report):
        fh.write('\t'.join(rxn_integration_report[cpd])+'\n')
fh.close()

with open(Biochem_Dir+'/'+Biochem+'_Mismatched_Compound_Integration_Report.txt','w') as fh:
    for cpd in sorted(mismatched_cpds, key=lambda k: len(mismatched_cpds[k]), reverse=True):
        fh.write('\t'.join([cpd,str(len(mismatched_cpds[cpd])),'|'.join(sorted(mismatched_cpds[cpd]))])+'\n')
fh.close()

print("Mismatched Compounds: "+str(len(mismatched_cpds.keys())))
#Here, for matches, re-write names, ecs and aliases
print("Saving additional ECs for "+str(len(New_EC_Count))+" reactions")
#reactions_helper.saveECs(ECs_Dict)
print("Saving additional names for "+str(len(New_Name_Count))+" reactions")
#reactions_helper.saveNames(Names_Dict)
print("Saving additional "+Biochem+" aliases for "+str(len(New_Alias_Count))+" reactions")
#reactions_helper.saveAliases(Original_Alias_Dict)
print("Saving "+str(len(New_Rxn_Count))+" new reactions from "+Biochem)
#reactions_helper.saveReactions(reactions_dict)

#Scripts to run afterwards
#./Rebuild_Reactions.py (in theory, because we adjust in this script, nothing should change)
#./Rebalance_Reactions.py (very important)
#./Adjust_Reaction_Protons.py
#./Adjust_Reaction_Water.py
#./Merge_Reactions.py (merges may happen because of water)
#./Update_Reaction_Aliases.py
