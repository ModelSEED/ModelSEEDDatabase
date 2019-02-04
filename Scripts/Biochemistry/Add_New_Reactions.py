#!/usr/bin/env python
import os, sys, re, copy
from csv import DictReader
from collections import OrderedDict
temp=list();
header=1;

Biochem="KEGG"
Biochem_Root="../../Biochemistry/Aliases/Provenance/Primary_Databases/";

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()
ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()
Reactions_Codes = ReactionsHelper.generateCodes(Reactions_Dict)

Default_Rxn = {"id":"cpd00001","name":"null","abbreviation":"null","aliases":"null",
               "code":"null","stoichiometry":"null","equation":"null","definition":"null",
               "reversibility":"=","direction":"=","deltag":"10000000","deltagerr":"10000000",
               "status":"NB","is_obsolete":0,"is_transport":0,
               "abstract_reaction":"null","pathways":"null","ec_numbers":"null",
               "compound_ids":"null","linked_reaction":"null","notes":"null","source":""}

Compounds_Alias_Dict=CompoundsHelper.loadMSAliases()
Source_Alias_Dict = dict()
for msid in Compounds_Alias_Dict.keys():
    for source in Compounds_Alias_Dict[msid].keys():
        if(source not in Source_Alias_Dict):
            Source_Alias_Dict[source]=dict()
        for alias in Compounds_Alias_Dict[msid][source]:
            if(alias not in Source_Alias_Dict[source]):
                Source_Alias_Dict[source][alias]=list()
            Source_Alias_Dict[source][alias].append(msid)

Original_Alias_Dict=ReactionsHelper.loadMSAliases()
New_Alias_Count=dict()

Names_Dict = ReactionsHelper.loadNames()
All_Names = dict()
New_Name_Count=dict()
for msid in sorted(Names_Dict.keys()):
    for name in Names_Dict[msid]:
        All_Names[name]=1

#Find last identifier and increment
last_identifier = list(sorted(Reactions_Dict.keys()))[-1]
identifier_count = int(re.sub('^rxn','',last_identifier))

New_Rxn_Count=dict()
Headers=list()
rxns=list()
missing_cpds=dict()
with open(Biochem_Root+Biochem+"_Reactions.tbl") as fh:
    for line in fh.readlines():
        line=line.strip()
        if(len(Headers)==0):
            Headers=line.split('\t')
            continue

        rxn=dict()
        array=line.split('\t',len(Headers))
        for i in range(len(Headers)):
            rxn[Headers[i]]=array[i]

        #Retrieve identifiers from within equation
        #Split based on whitespace, and remove compartment index
        original_cpd_array=rxn['EQUATION'].split(' ')
        new_cpd_array=list()
        for i in range(len(original_cpd_array)):
            if(re.search('\[[01]\]$',original_cpd_array[i])):
                new_cpd_array.append(re.sub('\[[01]\]$','',original_cpd_array[i]))

        all_matched=1
        for new_cpd in new_cpd_array:
            if(new_cpd in Source_Alias_Dict[Biochem]):
                msid = sorted(Source_Alias_Dict[Biochem][new_cpd])[0]
                rxn['EQUATION'] = re.sub(new_cpd,msid,rxn['EQUATION'])
            else:
                missing_cpds[new_cpd]=1
                all_matched=0

        if(all_matched==0):
            print("Warning: missing "+Biochem+" identifiers for reaction "+rxn['ID']+": "+rxn['EQUATION'])
            continue

        rxn_cpds_array = ReactionsHelper.parseEquation(rxn['EQUATION'])
        rxn_code = ReactionsHelper.generateCode(rxn_cpds_array)

        matched_rxn=None
        if(rxn_code in Reactions_Codes):
            matched_rxn = sorted(list(Reactions_Codes[rxn_code].keys()))[0]
            
            #Add Names and Alias
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

            #Add ID to aliases if the match is with a different reaction
            if(matched_rxn not in Original_Alias_Dict):
                Original_Alias_Dict[matched_rxn]=dict()
            if(matched_rxn in Original_Alias_Dict and Biochem not in Original_Alias_Dict[matched_rxn]):
                Original_Alias_Dict[matched_rxn][Biochem]=list()
            if(rxn['ID'] not in Original_Alias_Dict[matched_rxn][Biochem]):
                Original_Alias_Dict[matched_rxn][Biochem].append(rxn['ID'])
                New_Alias_Count[matched_rxn]=1

            #Update source type
            Reactions_Dict[matched_rxn]['source']='Primary Database'

        else:
            
            #New Reaction!
            #Generate new identifier
            identifier_count+=1
            new_identifier = 'rxn'+str(identifier_count)

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

            #Add source type
            new_rxn['source']='Primary Database'

            Reactions_Dict[new_rxn['id']]=new_rxn
            New_Rxn_Count[new_rxn['id']]=1

            #Rebuild key fields for reaction using parsed equation
            stoichiometry=ReactionsHelper.buildStoich(rxn_cpds_array)
            ReactionsHelper.rebuildReaction(Reactions_Dict[new_rxn['id']],stoichiometry)

            #Finally, because several new reactions may share equations
            if(rxn_code not in Reactions_Codes):
                Reactions_Codes[rxn_code]=dict()
            Reactions_Codes[rxn_code][new_rxn['id']]=1

print("Missing Compounds: "+"|".join(sorted(missing_cpds.keys())))

#Here, for matches, re-write names and aliases
print("Saving additional names for "+str(len(New_Name_Count))+" reactions")
ReactionsHelper.saveNames(Names_Dict)
print("Saving additional "+Biochem+" aliases for "+str(len(New_Alias_Count))+" reactions")
ReactionsHelper.saveAliases(Original_Alias_Dict)
print("Saving "+str(len(New_Rxn_Count.keys()))+" new reactions from "+Biochem)
ReactionsHelper.saveReactions(Reactions_Dict)

#Scripts to run afterwards
#./Merge_Reactions.py
#./Update_Reaction_Aliases.py
#./Rebalance_Reactions.py
#./Adjust_Reaction_Protons.py
#./Adjust_Reaction_Water.py
