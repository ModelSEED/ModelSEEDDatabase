#!/usr/bin/env python
import os, sys


columns={'2010':{'cpd_id':7,'formula':6,'structure':16,'rxn_id':9,'equation':8,'status':16},
         '2014':{'cpd_id':0,'formula':3,'structure':6,'rxn_id':0,'equation':6,'status':17}}

for year in ['2010','2014']:

    print("Statistics from "+year)
    compounds_dict=dict()
    with open(year+'/compoundTable.txt') as fh:
        for line in fh.readlines():
            line=line.strip()
            array=line.split('\t')
            compounds_dict[array[columns[year]['cpd_id']]]={'formula':array[columns[year]['formula']],
                                                            'structure':None}
            if(len(array)>16 and array[columns[year]['structure']].strip() != "" and array[columns[year]['structure']] != 'null'):
                compounds_dict[array[columns[year]['cpd_id']]]['structure']=array[columns[year]['structure']]
    fh.close()

    reactions_dict=dict()
    with open(year+'/reactionTable.txt') as fh:
        for line in fh.readlines():
            line=line.strip()
            array=line.split('\t')
            reactions_dict[array[columns[year]['rxn_id']]]={'eqn':array[columns[year]['equation']],'stat':array[columns[year]['status']]}
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
        print(entry,compound_counts[entry],float(compound_counts[entry])/float(compound_counts['cpd']))

    reaction_counts={'rxn':0,'Generic':0,'Complete':0,'Balanced':0}
    for rxn in reactions_dict:
        reaction_counts['rxn']+=1

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
        print(entry,reaction_counts[entry],float(reaction_counts[entry])/float(reaction_counts['rxn']))
    print("\n========================\n")
