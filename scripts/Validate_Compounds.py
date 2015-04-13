#! /usr/bin/python

import argparse
import re
from BiochemHelper import BiochemHelper

desc1 = '''
NAME
      Validate_Compounds -- validate a compounds file

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Validate a compounds file by checking for duplicate IDs and names, missing
      formulas, and invalid charges.  The cpdfile argument specifies the path
      to a compounds file which describes the compounds in a biochemistry.  The
      first line of the file is a header with the field names.  There is one
      compound per line with fields separated by tabs.

      The --charge optional argument specifies a value for invalid charges.  When
      the absolute value of the charge from the compound is larger than the value,
      the charge is invalid.  The default value is 50.

      When the --show-details optional argument is specified, details on all
      problems are displayed.  When the other --show optional arguments are
      specified, details on the corresponding type of problem are displayed.

      When the --fix-dup-names optional argument is specified, duplicate names
      are fixed by appending the string " (dupN)" to the name and abbreviation
      of duplicate compounds (where N is a number starting with 2).  The compounds
      file is rewritten with the fixed data.
'''

desc3 = '''
EXAMPLES
      Show summary data about a compounds file:
      > Validate_Compounds.py compounds.tsv
      
      Show details on compounds that have problems:
      > Validate_Compounds.py --show-details compounds.tsv

SEE ALSO
      Validate_Reactions

AUTHORS
      Mike Mundy 
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Validate_Compounds', epilog=desc3)
    parser.add_argument('cpdfile', help='path to compounds file', action='store')
    parser.add_argument('--charge', help='flag compounds with charge larger than value', action='store', dest='charge', type=int, default=50)
    parser.add_argument('--show-details', help='show details on all problems', action='store_true', dest='showDetails', default=False)
    parser.add_argument('--show-dup-ids', help='show details on duplicate IDs', action='store_true', dest='showDupIds', default=False)
    parser.add_argument('--show-bad-ids', help='show details on bad IDs', action='store_true', dest='showBadIds', default=False)
    parser.add_argument('--show-dup-names', help='show details on duplicate names', action='store_true', dest='showDupNames', default=False)
    parser.add_argument('--show-bad-names', help='show details on bad names', action='store_true', dest='showBadNames', default=False)
    parser.add_argument('--show-dup-abbrs', help='show details on duplicate abbreviations', action='store_true', dest='showDupAbbrs', default=False)
    parser.add_argument('--show-bad-abbrs', help='show details on bad abbreviations', action='store_true', dest='showBadAbbrs', default=False)
    parser.add_argument('--show-formulas', help='show details on missing formulas', action='store_true', dest='showFormulas', default=False)
    parser.add_argument('--show-charges', help='show details on invalid charges', action='store_true', dest='showCharges', default=False)
    parser.add_argument('--show-cofactors', help='show details on invalid cofactors', action='store_true', dest='showCofactors', default=False)
    parser.add_argument('--fix-dup-names', help='fix on duplicate names', action='store_true', dest='fixDupNames', default=False)
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
        args.showDupAbbrs = True
        args.showBadAbbrs = True
        args.showFormulas = True
        args.showCharges = True
        args.showCofactors = True

    # Read the compounds from the specified file.
    print 'Compound file: %s' %(args.cpdfile)
    helper = BiochemHelper()
    compounds = helper.readCompoundsFile(args.cpdfile)
    if compounds is None:
        print 'Error reading compounds file'
        exit(1)
    print 'Number of compounds: %d' %(len(compounds))

    # Check for duplicates, missing and invalid values.
    idDict = dict()
    duplicateId = 0
    badIdChars = list()
    nameDict = dict()
    duplicateName = 0
    badNameChars = list()
    abbrDict = dict()
    duplicateAbbr = 0
    badAbbrChars = list()
    noFormula = list()
    largeCharge = list()
    badCofactor = list()
    numCofactors = 0
    for index in range(len(compounds)):
        cpd = compounds[index]
        
        # Check for duplicate IDs.
        if cpd['id'] in idDict:
            if len(idDict[cpd['id']]) == 1:
                duplicateId += 1
            idDict[cpd['id']].append(index)
        else:
            idDict[cpd['id']] = [ index ]

        # Check for invalid characters in the ID.
        match = re.search(r'cpd\d\d\d\d\d', cpd['id'])
        if match is None:
            badIdChars.append(index)

        # Check for duplicate names.
        if cpd['name'] in nameDict:
            if len(nameDict[cpd['name']]) == 1:
                duplicateName += 1
            nameDict[cpd['name']].append(index)
        else:
            nameDict[cpd['name']] = [ index ]

        # Check for invalid characters in the name.
        try:
            cpd['name'].decode('ascii')
        except UnicodeDecodeError:
            badNameChars.append(index)

        # Check for duplicate abbreviations.
        if cpd['abbreviation'] in abbrDict:
            if len(abbrDict[cpd['abbreviation']]) == 1:
                duplicateAbbr += 1
            abbrDict[cpd['abbreviation']].append(index)
        else:
            abbrDict[cpd['abbreviation']] = [ index ]

        # Check for invalid characters in the abbreviation.
        try:
            cpd['abbreviation'].decode('ascii')
        except UnicodeDecodeError:
            badAbbrChars.append(index)

        # Check for missing or unknown formulas.
        if cpd['formula'] == '' or cpd['formula'] == 'noformula' or cpd['formula'] == 'unknown':
            noFormula.append(index)

        # Check for charges that are too big.
        if abs(cpd['defaultCharge']) > args.charge:
            largeCharge.append(index)

        # Check for invalid isCofactor flags.
        if cpd['isCofactor'] != 0 and cpd['isCofactor'] != 1:
            badCofactor.append(index)
        if cpd['isCofactor'] == 1:
            numCofactors += 1

    # Print summary data.
    print 'Number of compounds with duplicate IDs: %d' %(duplicateId)    
    print 'Number of compounds with bad characters in ID: %d' %(len(badIdChars))
    print 'Number of compounds with duplicate names: %d' %(duplicateName)
    print 'Number of compounds with bad characters in name: %d' %(len(badNameChars))
    print 'Number of compounds with duplicate abbreviations: %d' %(duplicateAbbr)
    print 'Number of compounds with bad characters in abbreviation: %d' %(len(badAbbrChars))
    print 'Number of compounds with no formula: %d' %(len(noFormula))
    print 'Number of compounds with charge larger than %d: %d' %(args.charge, len(largeCharge))
    print 'Number of compounds with bad isCofactor flag: %d' %(len(badCofactor))
    print 'Number of compounds flagged as cofactor: %d' %(numCofactors)
    print

    # Print details if requested.
    if args.showDupIds:
        for id in idDict:
            if len(idDict[id]) > 1:
                print 'Duplicate compound ID: %s' %(id)
                for dup in idDict[id]:
                    print 'Line %05d: %s' %(compounds[dup]['linenum'], compounds[dup])
                print
    if args.showBadIds:
        if len(badIdChars) > 0:
            print 'Compounds with bad characters in ID:'
            for index in range(len(badIdChars)):
                print 'Line %05d: %s' %(compounds[badIdChars[index]]['linenum'], compounds[badIdChars[index]])
            print
    if args.showDupNames:
        for name in nameDict:
            if len(nameDict[name]) > 1:
                print 'Duplicate compound name: %s' %(name)
                for dup in nameDict[name]:
                    print 'Line %05d: %s' %(compounds[dup]['linenum'], compounds[dup])
                print
    if args.showBadNames:
        if len(badNameChars) > 0:
            print 'Compounds with bad characters in name:'
            for index in range(len(badNameChars)):
                print 'Line %05d: %s' %(compounds[badNameChars[index]]['linenum'], compounds[badNameChars[index]])
            print
    if args.showDupAbbrs:
        for abbr in abbrDict:
            if len(abbrDict[abbr]) > 1:
                print 'Duplicate compound abbreviation: %s' %(abbr)
                for dup in abbrDict[abbr]:
                    print 'Line %05d: %s' %(compounds[dup]['linenum'], compounds[dup])
                print
    if args.showBadAbbrs:
        if len(badAbbrChars) > 0:
            print 'Compounds with bad characters in abbreviation:'
            for index in range(len(badAbbrChars)):
                print 'Line %05d: %s' %(compounds[badAbbrChars[index]]['linenum'], compounds[badAbbrChars[index]])
            print
    if args.showFormulas:
        if len(noFormula) > 0:
            print 'Compounds with no formula:'
            for index in range(len(noFormula)):
                print 'Line %05d: %s' %(compounds[noFormula[index]]['linenum'], compounds[noFormula[index]])
            print
    if args.showCharges:
        if len(largeCharge) > 0:
            print 'Compounds with charge larger than %d:' %(args.charge)
            for index in range(len(largeCharge)):
                print 'Line %05d: %s' %(compounds[largeCharge[index]]['linenum'], compounds[largeCharge[index]])
            print
    if args.showCofactors:
        if len(badCofactor) > 0:
            print 'Compounds with bad isCofactor flag:'
            for index in range(len(badCofactor)):
                print 'Line %05d: %s' %(compounds[badCofactor[index]]['linenum'], compounds[badCofactor[index]])
            print

    if args.fixDupNames:
        for name in nameDict:
            if len(nameDict[name]) > 1:
                for index in range(1,len(nameDict[name])): # Leave the first duplicate unchanged
                    dup = nameDict[name][index]
                    compounds[dup]['name'] += ' (dup%d)' %(index+1)
                    compounds[dup]['abbreviation'] = compounds[dup]['name']
                    if compounds[dup]['formula'] != compounds[nameDict[name][0]]['formula']:
                        print 'WARNING: formula mismatch'
                        print 'Line %05d: %s' %(compounds[nameDict[name][0]]['linenum'], compounds[nameDict[name][0]])
                        print 'Line %05d: %s' %(compounds[dup]['linenum'], compounds[dup])

        with open(args.cpdfile, 'w') as handle:
            handle.write('id\tname\tabbreviation\tformula\tcharge\tisCofactor\n')
            for index in range(len(compounds)):
                cpd = compounds[index]
                line = '%s\t%s\t%s\t%s\t%d\t%d\n' %(cpd['id'], cpd['name'], cpd['abbreviation'], cpd['formula'], cpd['defaultCharge'], cpd['isCofactor'])
                handle.write(line)

    exit(0)
