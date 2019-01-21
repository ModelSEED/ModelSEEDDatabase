import sys
import os
import json
import datetime
import argparse
from Bio.ExPASy import Enzyme
import xmltodict
import pprint

timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description = '', formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--test', '-t',
    action = 'store_true',
    help = 'Use local flat files rather than pulling from external source' )

parser.add_argument('--summary', '-s',
    action = 'store_true',
    help = 'Print summary info' )

args = parser.parse_args()

## FUNCTIONS ###################################################################

def get_KEGG_KOs():

    ko_count = 0
    ko_ontologyDictionary_filename = 'KEGG_KO_ontologyDictionary.json'

    # get KEGG version number
    kegg_release_call = 'curl -g -s -S http://rest.kegg.jp/info/kegg'
    kegg_release_raw = os.popen(kegg_release_call).read()
    kegg_release = kegg_release_raw.split("\n")[1].split()[2].rstrip(",")

    # create dictionary
    kegg_dict = {'data_version'  : kegg_release,
                'date'           : timestamp,
                'format_version' : 'N/A',
                'ontology'       : 'kegg_orthology',
                'term_hash'      : {}
                }

    # get raw orthology from API or local file
    kegg_raw = ""

    if args.test == False:
        kegg_api_call = 'curl -g -s -S http://rest.kegg.jp/list/orthology'
        kegg_raw = os.popen(kegg_api_call).read()
    else:
        with open('kegg_orthologs.txt', 'r') as myfile:
            kegg_raw = myfile.read()

    # parse orthology data
    for line in kegg_raw.split("\n"):

        if len(line) > 0:
            id, name_raw = line.split("\t")
            id = id.replace("ko:", "")

            nameSplit = name_raw.split(";")

            name = ""
            synonym = ""

            if len(nameSplit) > 1:
                name = nameSplit[1].strip()
                synonym = nameSplit[0].strip()
            else:
                print(name)

            kegg_dict['term_hash'][id] = {'id' : id, 'name' : name, 'synonyms' : synonym.split(", ")}
            ko_count += 1

    # save json file
    with open(ko_ontologyDictionary_filename, 'w') as outfile:
        json.dump(kegg_dict, outfile, indent = 2)

    # print summary
    if args.summary == True:
        print("kegg_orthology", ko_count, kegg_release, ko_ontologyDictionary_filename, sep = "\t")

def get_KEGG_RXNs():

    rxn_count = 0
    rxn_ontologyDictionary_filename = 'KEGG_RXN_ontologyDictionary.json'

    # get KEGG version number
    kegg_release_call = 'curl -g -s -S http://rest.kegg.jp/info/kegg'
    kegg_release_raw = os.popen(kegg_release_call).read()
    kegg_release = kegg_release_raw.split("\n")[1].split()[2].rstrip(",")

    # create dictionary
    kegg_dict = {'data_version'  : kegg_release,
                'date'           : timestamp,
                'format_version' : 'N/A',
                'ontology'       : 'kegg_reactions',
                'term_hash'      : {}
                }

    # get raw orthology from API or local file
    kegg_raw = ""

    if args.test == False:
        kegg_api_call = 'curl -g -s -S http://rest.kegg.jp/list/reaction'
        kegg_raw = os.popen(kegg_api_call).read()
    else:
        with open('kegg_reactions.txt', 'r') as myfile:
            kegg_raw = myfile.read()

    # parse reaction data
    for line in kegg_raw.split("\n"):

        if len(line) > 0:

            id, name_raw = line.split("\t")
            id = id.replace("rn:", "")

            fullLineSplit = name_raw.split("; ")

            name = fullLineSplit.pop(0)

            kegg_dict['term_hash'][id] = {'id' : id, 'name' : name, 'synonyms' : fullLineSplit}
            rxn_count += 1

    # save json file
    with open(rxn_ontologyDictionary_filename, 'w') as outfile:
        json.dump(kegg_dict, outfile, indent = 2)

    # print summary
    if args.summary == True:
        print("kegg_reactions", rxn_count, kegg_release, rxn_ontologyDictionary_filename, sep = "\t")

def get_EC_RXNs():

    ec_count = 0
    ec_ontologyDictionary_filename = 'EBI_EC_ontologyDictionary.json'

    # get EC version number
    ebi_ec_release = ""

    if args.test == False:
        ebi_ec_call = 'curl ftp://ftp.ebi.ac.uk/pub/databases/enzyme/enzclass.txt'
        ebi_ec_release = os.popen(ebi_ec_call).read()
    else:
        with open('enzclass.txt', 'r') as myfile:
            ebi_ec_release = myfile.read()

    ebi_ec_release = ebi_ec_release.split("\n")[7].split()[1]

    # create dictionary
    ec_dict = {'data_version'    : ebi_ec_release,
                'date'           : timestamp,
                'format_version' : 'N/A',
                'ontology'       : 'ec_orthology',
                'term_hash'      : {}
                }

    # parse data
    ebi_ec_enzyme = 'enzyme.dat'

    if args.test == False:
        ebi_ec_call = 'curl ftp://ftp.ebi.ac.uk/pub/databases/enzyme/enzyme.dat > enzyme.dat'
        os.system(ebi_ec_call)

    records = Enzyme.parse(open(ebi_ec_enzyme))

    for record in records:
        ec_dict['term_hash'][record['ID']] = {'id'       : record['ID'],
                                              'name'     : record['DE'],
                                              'synonyms' : record['AN']
                                             }
        ec_count += 1

    # save json file
    with open(ec_ontologyDictionary_filename, 'w') as outfile:
        json.dump(ec_dict, outfile, indent = 2)

    # print summary
    if args.summary == True:
        print("ec_orthology", ec_count, ebi_ec_release, ec_ontologyDictionary_filename, sep = "\t")

def get_METACYC_RXNs():

    # get_METACYC_RXNs() does not have an option to pull from an online database, only from the local MetaCycReactionsFull.xml file

    metacyc_count = 0
    metacyc_ontologyDictionary_filename = 'METACYC_ontologyDictionary.json'

    # convert to json and upload
    with open('MetaCycReactionsFull.xml') as fd:
        doc = xmltodict.parse(fd.read())
    with open('metacyc_temp.json', 'w') as outfile:
        json.dump(doc, outfile, indent = 2)
    metacyc_raw = json.loads(open("metacyc_temp.json", "r").read() )

    # get metacyc version number
    metacyc_version = metacyc_raw['ptools-xml']['@ptools-version']

    # create dictionary
    metacyc_dict = {'data_version'    : metacyc_version,
                    'date'            : timestamp,
                    'format_version'  : 'N/A',
                    'ontology'        : 'metacyc_orthology',
                    'term_hash'       : {}
                   }

    # parse file
    for reaction in metacyc_raw['ptools-xml']['Reaction']:

        id = reaction['@ID']
        nameList = []

        if 'ec-number' in reaction: # also removes the "EC-"
            if type(reaction['ec-number']) == dict:
                nameList.append(reaction['ec-number']['#text'].replace("EC-", ""))
            elif type(reaction['ec-number']) == str:
                nameList.append(reaction['ec-number'].replace("EC-", ""))

        if 'enzymatic-reaction' in reaction:

            for level1 in reaction['enzymatic-reaction']['Enzymatic-Reaction']:

                if type(level1) == dict:

                    # some, like META:1TRANSKETO-RXN, don't have an @ID everywhere... so they'll get skipped
                    if '@ID' in level1:
                        nameList.append(level1['@ID'])

                        for subreaction in level1:
                            if 'common-name' in level1:
                                nameList.append(level1['common-name']['#text'])

                elif type(level1) == str:

                    # some, like THYROID-HORMONE-AMINOTRANSFERASE-RXN, don't seem to have a common-name, so they'll get skipped
                    if 'common-name' in reaction['enzymatic-reaction']['Enzymatic-Reaction']:
                        nameList.append(reaction['enzymatic-reaction']['Enzymatic-Reaction']['common-name']['#text'])

        nameList = list(set(nameList))

        # fix unicode, but doesn't currently work
        cleanNameList = []
        for name in nameList:
            #name = name.replace("&beta;", "β")
            cleanNameList.append(name)

        # save the first item from the name list as the name, the rest are synonyms
        if len(cleanNameList) > 0:
            name = cleanNameList.pop(0)
        else:
            name = "NONE"

        metacyc_dict['term_hash'][id] = {'id'       : id,
                                         'name'     : name,
                                         'synonyms' : cleanNameList
                                         }

        metacyc_count += 1

    # save json file
    with open(metacyc_ontologyDictionary_filename, 'w') as outfile:
        json.dump(metacyc_dict, outfile, indent = 2)

    # print summary
    if args.summary == True:
        print("metacyc_orthology", metacyc_count, metacyc_version, metacyc_ontologyDictionary_filename, sep = "\t")

## RUN FUNCTIONS ###############################################################

get_KEGG_KOs()
get_KEGG_RXNs()
get_EC_RXNs()
get_METACYC_RXNs()
