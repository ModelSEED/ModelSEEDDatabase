#! /usr/bin/python

import argparse
import re
from BiochemHelper import BiochemHelper

desc1 = '''
NAME
      Validate_Reactions -- validate a reactions file

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Validate a reactions file by checking for duplicate IDs and names, missing
      equations, and invalid status.  The rxnfile argument specifies the path
      to a reactions file which describes the reactions in a biochemistry.  The
      first line of the file is a header with the field names.  There is one
      reaction per line with fields separated by tabs.

      When the --show-details optional argument is specified, details on all
      problems are displayed.  When the other --show optional arguments are
      specified, details on the corresponding type of problem are displayed.
'''

desc3 = '''
EXAMPLES
      Show summary data about a reactions file:
      > Validate_Reactions.py reactions.tsv
      
      Show details on reactions that have problems:
      > Validate_Reactions.py --show-details reactions.tsv

SEE ALSO
      Validate_Compounds

AUTHORS
      Mike Mundy 
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Validate_Reactions', epilog=desc3)
    parser.add_argument('rxnfile', help='path to reactions file', action='store')
    parser.add_argument('--show-details', help='show details on all problems', action='store_true', dest='showDetails', default=False)
    parser.add_argument('--show-dup-ids', help='show details on duplicate IDs', action='store_true', dest='showDupIds', default=False)
    parser.add_argument('--show-bad-ids', help='show details on bad IDs', action='store_true', dest='showBadIds', default=False)
    parser.add_argument('--show-dup-names', help='show details on duplicate names', action='store_true', dest='showDupNames', default=False)
    parser.add_argument('--show-bad-names', help='show details on bad names', action='store_true', dest='showBadNames', default=False)
    parser.add_argument('--show-bad-eq', help='show details on missing reactants or products in equations', action='store_true', dest='showBadEquations', default=False)
    parser.add_argument('--show-status', help='show details on status types', action='store_true', dest='showStatus', default=False)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # The --show-details option turns on all of the detail options.
    if args.showDetails:
        args.showDupIds = True
        args.showBadIds = True
        args.showDupNames = True
        args.showBadNames = True
        args.showStatus = True

    # Read the reactions from the specified file.
    helper = BiochemHelper()
    reactions = helper.readReactionsFile(args.rxnfile)
    print 'Reaction file: %s' %(args.rxnfile)
    print 'Number of reactions: %d' %(len(reactions))

    # Check for duplicates, missing and invalid values.
    idDict = dict()
    duplicateId = 0
    badIdChars = list()
    nameDict = dict()
    duplicateName = 0
    badNameChars = list()
    noEquation = list()
    noDefinition = list()
    noReactants = list()
    noProducts = list()
    statusTypes = dict()
    okStatus = 0
    for index in range(len(reactions)):
        rxn = reactions[index]
        
        # Check for duplicate IDs.
        if rxn['id'] in idDict:
            if len(idDict[rxn['id']]) == 1:
                duplicateId += 1
            idDict[rxn['id']].append(index)
        else:
            idDict[rxn['id']] = [ index ]

        # Check for invalid characters in the ID.
        match = re.search(r'rxn\d\d\d\d\d', rxn['id'])
        if match is None:
            badIdChars.append(index)

        # Check for duplicate names.
        if rxn['name'] in nameDict:
            if len(nameDict[rxn['name']]) == 1:
                duplicateName += 1
            nameDict[rxn['name']].append(index)
        else:
            nameDict[rxn['name']] = [ index ]

        # Check for invalid characters in the name.
        try:
            rxn['name'].decode('ascii')
        except UnicodeDecodeError:
            badNameChars.append(index)

        # Check for missing reactants and/or products.
        reactants, products = helper.parseEquation(rxn['equation'])
        if reactants is None and products is None:
            noEquation.append(index)
        if len(reactants) == 0:
            noReactants.append(index)
        if len(products) == 0:
            noProducts.append(index)

        reactants, products = helper.parseEquation(rxn['definition'])
        if reactants is None and products is None:
            noDefinition.append(index)

        # Check reaction status.
        if rxn['status'] == 'OK':
            okStatus += 1
        else:
            fields = rxn['status'].split('|')
            for index in range(len(fields)):
                if ':' in fields[index]:
                    pos = fields[index].find(':')
                    type = fields[index][:pos]
                else:
                    type = fields[index]
                if type in statusTypes:
                    statusTypes[type] += 1
                else:
                    statusTypes[type] = 1

    # Print summary data.
    print 'Number of reactions with duplicate IDs: %d' %(duplicateId)    
    print 'Number of reactions with bad characters in ID: %d' %(len(badIdChars))
    print 'Number of reactions with duplicate names: %d' %(duplicateName)
    print 'Number of reactions with bad characters in name: %d' %(len(badNameChars))
    print 'Number of reactions with missing equation: %d' %(len(noEquation))
    print 'Number of reactions with missing definition: %d' %(len(noDefinition))
    print 'Number of reactions with no reactants: %d' %(len(noReactants))
    print 'Number of reactions with no products: %d' %(len(noProducts))
    print 'Number of reactions with OK status: %d' %(okStatus)
    print 'Number of reactions with error status: %d' %(len(reactions)-okStatus)
    print

    # Print details if requested.
    if args.showDupIds:
        for id in idDict:
            if len(idDict[id]) > 1:
                print 'Duplicate reaction ID: %s' %(id)
                for dup in idDict[id]:
                    print 'Line %05d: %s' %(reactions[dup]['lineno'], reactions[dup])
                print
    if args.showBadIds:
        if len(badIdChars) > 0:
            print 'reactions with bad characters in ID:'
            for index in range(len(badIdChars)):
                print 'Line %05d: %s' %(reactions[badIdChars[index]]['lineno'], reactions[badIdChars[index]])
            print
    if args.showDupNames:
        for name in nameDict:
            if len(nameDict[name]) > 1:
                print 'Duplicate reaction name: %s' %(name)
                for dup in nameDict[name]:
                    print 'Line %05d: %s' %(reactions[dup]['lineno'], reactions[dup])
                print
    if args.showBadNames:
        if len(badNameChars) > 0:
            print 'reactions with bad characters in name:'
            for index in range(len(badNameChars)):
                print 'Line %05d: %s' %(reactions[badNameChars[index]]['lineno'], reactions[badNameChars[index]])
            print
    if args.showBadEquations:
        if len(noReactants) > 0:
            print 'Reactions with no reactants:'
            for index in range(len(noReactants)):
                print 'Line %05d: %s' %(reactions[noReactants[index]]['lineno'], reactions[noReactants[index]])
            print
        if len(noProducts) > 0:
            print 'Reactions with no products:'
            for index in range(len(noProducts)):
                print 'Line %05d: %s' %(reactions[noProducts[index]]['lineno'], reactions[noProducts[index]])
            print
        if len(noEquation) > 0:
            print 'Reactions with no equation:'
            for index in range(len(noEquation)):
                print 'Line %05d: %s' %(reactions[noEquation[index]]['lineno'], reactions[noEquation[index]])
            print
        if len(noDefinition) > 0:
            print 'Reactions with no definition:'
            for index in range(len(noDefinition)):
                print 'Line %05d: %s' %(reactions[noDefinition[index]]['lineno'], reactions[noDefinition[index]])
            print
    if args.showStatus:
        for type in statusTypes:
            print 'Reaction status %s: %d' %(type, statusTypes[type])

    exit(0)
