#! /usr/bin/python

import argparse
import re

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

      When the --show-details optional argument is specified, details on each
      duplicate or missing value are displayed.
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
    parser.add_argument('--show-details', help='show details on problems', action='store_true', dest='showDetails', default=False)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Read the reactions from the specified file.
    reactions = list()
    with open(args.rxnfile, 'r') as handle:
        # The first line has the header with the field names.
        fieldNames = handle.readline().strip().split('\t')
        # @todo Validate required fields are in the header
        
        lineno = 1
        for line in handle:
            lineno += 1
            fields = line.strip().split('\t')
            if len(fields) < len(fieldNames):
                print 'WARNING: Reaction on line %d is missing one or more fields, %s' %(lineno, fields)
                continue
            rxn = dict()
            rxn['id'] = fields[0]
            rxn['name'] = fields[1]
            rxn['equation'] = fields[2]
            rxn['definition'] = fields[3]
            rxn['status'] = fields[4]
            rxn['lineno'] = lineno
            reactions.append(rxn)

    print 'Reaction file: %s' %(args.rxnfile)
    print 'Number of reactions: %d' %(len(reactions))

    # Check for duplicates, missing and invalid values.
    idDict = dict()
    duplicateId = 0
    badIdChars = list()
    nameDict = dict()
    duplicateName = 0
    badNameChars = list()
    noEquation = 0
    noDefinition = 0
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

        # Check for missing equation or definition.
        if len(rxn['equation']) == 0:
            noEquation += 1
        if len(rxn['definition']) == 0:
            noDefinition += 1

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
    print 'Number of reactions with missing equation: %d' %(noEquation)
    print 'Number of reactions with missing definition: %d' %(noDefinition)
    print 'Number of reactions with OK status: %d' %(okStatus)
    print 'Number of reactions with error status: %d' %(len(reactions)-okStatus)
    print

    # Print details if requested.
    if args.showDetails:
        for id in idDict:
            if len(idDict[id]) > 1:
                print 'Duplicate reaction ID: %s' %(id)
                for dup in idDict[id]:
                    print 'Line %05d: %s' %(reactions[dup]['lineno'], reactions[dup])
                print
        if len(badIdChars) > 0:
            print 'reactions with bad characters in ID:'
            for index in range(len(badIdChars)):
                print 'Line %05d: %s' %(reactions[badIdChars[index]]['lineno'], reactions[badIdChars[index]])
            print
        for name in nameDict:
            if len(nameDict[name]) > 1:
                print 'Duplicate reaction name: %s' %(name)
                for dup in nameDict[name]:
                    print 'Line %05d: %s' %(reactions[dup]['lineno'], reactions[dup])
                print
        if len(badNameChars) > 0:
            print 'reactions with bad characters in name:'
            for index in range(len(badNameChars)):
                print 'Line %05d: %s' %(reactions[badNameChars[index]]['lineno'], reactions[badNameChars[index]])
            print
        for type in statusTypes:
            print 'Reaction status %s: %d' %(type, statusTypes[type])

    exit(0)
