#!/usr/bin/env python
import os, sys


columns={'2010':{'cpd_id':7,'formula':6,'structure':16,'rxn_id':9,'equation':8,'status':16,'reversibility':18,'direction':14},
         '2014':{'cpd_id':0,'formula':3,'structure':6,'rxn_id':0,'equation':6,'status':17,'reversibility':8,'direction':9}}

out = open('Growth_Stats.tsv','w')
Old_Rxn_Rev=dict()
for year in ['2010','2014']:

    print("Statistics from "+year)
    compounds_dict=dict()
    with open('../'+year+'/compoundTable.txt') as fh:
        for line in fh.readlines():
            line=line.strip()
            array=line.split('\t')
            compounds_dict[array[columns[year]['cpd_id']]]={'formula':array[columns[year]['formula']],
                                                            'structure':None}
            if(len(array)>16 and array[columns[year]['structure']].strip() != "" and array[columns[year]['structure']] != 'null'):
                compounds_dict[array[columns[year]['cpd_id']]]['structure']=array[columns[year]['structure']]
    fh.close()

    reactions_dict=dict()
    with open('../'+year+'/reactionTable.txt') as fh:
        for line in fh.readlines():
            line=line.strip()
            array=line.split('\t')
            reactions_dict[array[columns[year]['rxn_id']]]={'eqn':array[columns[year]['equation']],
                                                            'stat':array[columns[year]['status']]}

            if(year=='2014'):
                Old_Rxn_Rev[array[columns[year]['rxn_id']]]=array[columns[year]['reversibility']]

    fh.close()

    compound_counts={'cpd':0,'Generic':0,'Structure':0}
    for cpd in compounds_dict:
        compound_counts['cpd']+=1
        if(compounds_dict[cpd]['structure'] is not None):
            compound_counts['Structure']+=1
        if('R' in compounds_dict[cpd]['formula']):
            compound_counts['Generic']+=1

    print(str(compound_counts['cpd'])+" compounds")
    for entry in ['Structure','Generic']:
        pct = "{0:.2f}".format(float(compound_counts[entry])/float(compound_counts['cpd']))
        print(entry,compound_counts[entry],pct)
    out.write('\t'.join([year,'cpd',str(compound_counts['cpd']),str(compound_counts['Structure'])])+'\n')

    reaction_counts={'rxn':0,'Generic':0,'Complete':0,'Balanced':0}
    for rxn in reactions_dict:
        reaction_counts['rxn']+=1

        #Skip empty reactions
        if(reactions_dict[rxn]['eqn'] == " <=> "):
            continue

        complete=True
        generic=False

        for entry in reactions_dict[rxn]['eqn'].split(' '):
            entry=entry.strip()

            if(len(entry)>0 and entry[-1] == ']'):
                entry=entry[:-3]

            if(len(entry)>0 and entry[0:3] == 'n-1'):
                entry=entry[3:]

            if(entry in compounds_dict):

                if(compounds_dict[entry]['structure'] is None):
                    complete=False

                if('R' in compounds_dict[entry]['formula']):
                    generic=True

        if(complete is True):
            reaction_counts['Complete']+=1

        if(generic is True):
            reaction_counts['Generic']+=1

        if('OK' in reactions_dict[rxn]['stat']):
            reaction_counts['Balanced']+=1

    print(str(reaction_counts['rxn'])+" reactions")
    for entry in ['Complete','Balanced','Generic']:
        pct = "{0:.2f}".format(float(reaction_counts[entry])/float(reaction_counts['rxn']))
        print(entry,reaction_counts[entry],pct)
    print("\n========================\n")
    out.write('\t'.join([year,'rxn',str(reaction_counts['rxn']),str(reaction_counts['Balanced'])])+'\n')

from BiochemPy import Compounds, Reactions,InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

print("Statistics from 2019")
out = open('Growth_Stats.tsv','a')
compound_counts={'Generic':0,'Structure':0}
for cpd in compounds_dict:
    if(compounds_dict[cpd]['smiles'] != ""):
        compound_counts['Structure']+=1
    if('R' in compounds_dict[cpd]['formula']):
        compound_counts['Generic']+=1

print(str(len(compounds_dict.keys()))+" compounds")
for entry in ['Structure','Generic']:
    pct = "{0:.2f}".format(float(compound_counts[entry])/float(len(compounds_dict.keys())))
    print(entry,compound_counts[entry],pct)
out.write('\t'.join(['2019','cpd',str(len(compounds_dict.keys())),str(compound_counts['Structure'])])+'\n')

reaction_counts={'Generic':0,'Complete':0,'Balanced':0}
New_Rxn_Rev=dict()
for rxn in reactions_dict:

    if(reactions_dict[rxn]['status'] == "EMPTY"):
        continue

    complete=True
    generic=False

    New_Rxn_Rev[rxn]=reactions_dict[rxn]['reversibility']

    for entry in reactions_dict[rxn]['compound_ids'].split(';'):

        if(entry in compounds_dict):

            if(compounds_dict[entry]['smiles'] == ""):
                complete=False

            if('R' in compounds_dict[entry]['formula']):
                generic=True

    if(complete is True):
        reaction_counts['Complete']+=1

    if(generic is True):
        reaction_counts['Generic']+=1

    if('OK' in reactions_dict[rxn]['status']):
        reaction_counts['Balanced']+=1

print(str(len(reactions_dict.keys()))+" reactions")
for entry in ['Complete','Balanced','Generic']:
    pct = "{0:.2f}".format(float(reaction_counts[entry])/float(len(reactions_dict.keys())))
    print(entry,reaction_counts[entry],pct)
out.write('\t'.join(['2019','rxn',str(len(reactions_dict.keys())),str(reaction_counts['Balanced'])])+'\n')
out.close()

print("\n========================")
print("For Thermodynamics section")
Rev_Counts={'New_Rev':0,'New_iRev':0,'ReviRev':0,'iRevRev':0}
for rxn in New_Rxn_Rev:
    if(New_Rxn_Rev[rxn] =='?'):
        continue

    if(rxn not in Old_Rxn_Rev):
        if(New_Rxn_Rev[rxn]=='='):
            Rev_Counts['New_Rev']+=1
        else:
            Rev_Counts['New_iRev']+=1
    else:
        if(Old_Rxn_Rev[rxn] == '?'):
            if(New_Rxn_Rev[rxn]=='='):
                Rev_Counts['New_Rev']+=1
            else:
                Rev_Counts['New_iRev']+=1
        else:
            if(Old_Rxn_Rev[rxn] == '=' and New_Rxn_Rev[rxn] != '='):
                Rev_Counts['ReviRev']+=1
            if(Old_Rxn_Rev[rxn] != '=' and New_Rxn_Rev[rxn] == '='):
                Rev_Counts['iRevRev']+=1

print("Reversibility for New Reactions: ",Rev_Counts['New_Rev']+Rev_Counts['New_iRev'])
print("Reversible to irreversible: ",Rev_Counts['ReviRev'])
print("Irreversible to reversible: ",Rev_Counts['iRevRev'])
