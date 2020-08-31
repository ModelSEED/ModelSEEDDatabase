#!/usr/bin/env python
import os, sys, re
temp=list();
from BiochemPy import Reactions, Compounds
compounds_helper=Compounds()
compounds_dict = compounds_helper.loadCompounds()
reactions_helper=Reactions()
reactions_dict = reactions_helper.loadReactions()
reactions_codes = reactions_helper.generateCodes(reactions_dict,stoich=False,transport=False)

names_dict = compounds_helper.loadNames()
searchnames_dict = dict()
for msid in sorted(names_dict):
    for name in names_dict[msid]:
        searchname = compounds_helper.searchname(name)
        #Avoid redundancy where possible
        if(searchname not in searchnames_dict):
            searchnames_dict[searchname]=msid

mhc = open('Mishit_Compound_Names.txt','w')
with open('Parsed_Enzyme_Equations.txt') as fh:
    for line in fh.readlines():
        line=line.strip()
        (id,old_equation)=line.split('\t')
        array = re.split(' (<?=>?|\+) ', old_equation)

        new_array=list()
        mishit=False
        for i in range(len(array)):
            if(array[i] == '+' or array[i] == '='):
                new_array.append(array[i])
                continue

            array[i] = array[i].strip()
            cpt="0"

            #Remove but preserving transport indicators
            match = re.search('\(Side\s?([12])?\)',array[i])
            if(match is not None):
                if(match.groups(0)[0]=='2'):
                    cpt="1"
                array[i] = re.sub('\(Side\s?[12]?\)','',array[i])
            match = re.search('\((In|Out)\)',array[i])
            if(match is not None):
                if(match.groups(0)[0]=='Out'):
                    cpt="1"
                array[i] = re.sub('\((In|Out)\)','',array[i])

            #Remove stoichiometry
            array[i] = re.sub('^\(?(\d?-)?[\dn]+\)? ','',array[i])
            array[i] = re.sub('^\d ','',array[i]) #Do it again because previous regex somehow misses single digits in parentheses

            #(n), (n+1), (n-1), (n-2),(n-x),(m), (n+m)
            array[i] = re.sub('\([nm]([+-][12xm])?\)','',array[i])

            #(n-1 to 5)
            array[i] = array[i].replace('(n-1 to 5)','')

            #(omega=180)
            array[i] = array[i].replace('(omega=180)','')

            #(65), (45-87), (27/28)
            array[i] = re.sub('\((\d{2,}[\-\/])?\d{2,}\)','',array[i])

            #want to remove all instances of '(\d)$' except for when its in a formula
            #H(2), AH(2), NH(3), N(2), O(2), CO(2), R'NH(2), RCH(2)NH(2), ArCH(2)NH(2), H(2)CO(3), H(2)O(2)
            match = re.search('(((^|[\'\)])[CAN\s]?[NHO]))\(\d\)$', array[i])
            if(match is None):
                array[i] = re.sub('\(\d\)$','',array[i])

            searchname = compounds_helper.searchname(array[i])
            if(searchname in searchnames_dict):
                new_array.append(searchnames_dict[searchname]+'['+cpt+']')
            else:
                mishit=True
                mhc.write(array[i]+'\n')

        if(mishit is False):
            equation_string = " ".join(new_array)
            rxn_cpds_array = reactions_helper.parseEquation(equation_string)
            new_rxn_cpds_array = reactions_helper.removeCpdRedundancy(rxn_cpds_array)
            rxn_code = reactions_helper.generateCode(new_rxn_cpds_array,stoich=False,transport=False)
            if(rxn_code in reactions_codes):
                for rxn in reactions_codes[rxn_code]:
                    print(rxn, id, old_equation, equation_string)
                    pass
