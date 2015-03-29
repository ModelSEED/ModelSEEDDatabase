#! /usr/bin/python

import argparse
from BiochemHelper import BiochemHelper

desc1 = '''
NAME
      Compare_Compounds -- compare two compounds files

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Compare two compounds files and find differences in compounds that have
      have the same ID.  The cpdfile1 and cpdfile2 arguments specify the paths
      to compounds files which describe the compounds in a biochemistry.  The
      first line of the file is a header with the field names.  There is one
      compound per line with fields separated by tabs.

      When the --show-names optional argument is specified, details on compounds
      with different names are displayed.  When the --show-formulas optional
      argument is specified, details on compounds with different formulas are
      displayed.  When the --show-charges optional argument is specified, details
      on compounds with different charges are displayed.  When the --show-details
      optional argument is specified, all differences are displayed.
'''

desc3 = '''
EXAMPLES
      Show summary data about differences between two compounds files:
      > Compare_Compounds.py compounds1.tsv compounds2.tsv
      
SEE ALSO
      Compare_Reactions

AUTHORS
      Mike Mundy 
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Compare_Compounds', epilog=desc3)
    parser.add_argument('cpdfile1', help='path to first compounds file', action='store')
    parser.add_argument('cpdfile2', help='path to second compounds file', action='store')
    parser.add_argument('--show-details', help='show details on all differences', action='store_true', dest='showDetails', default=False)
    parser.add_argument('--show-names', help='show details on compounds with different names', action='store_true', dest='showNames', default=False)
    parser.add_argument('--show-formulas', help='show details on compounds with different formulas', action='store_true', dest='showFormulas', default=False)
    parser.add_argument('--show-charges', help='show details on compounds with different charges', action='store_true', dest='showCharges', default=False)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # The --show-details option turns on all of the detail options.
    if args.showDetails:
        args.showNames = True
        args.showFormulas = True
        args.showCharges = True

    # Read the compounds from the first file.
    helper = BiochemHelper()
    firstCompounds = helper.buildDictFromListOfObjects(helper.readCompoundsFile(args.cpdfile1))
    print 'First compound file: %s' %(args.cpdfile1)
    print '  Number of compounds: %d' %(len(firstCompounds))

    # Read the compounds from the second file.
    secondCompounds = helper.buildDictFromListOfObjects(helper.readCompoundsFile(args.cpdfile2))
    print 'Second compound file: %s' %(args.cpdfile2)
    print '  Number of compounds: %d' %(len(secondCompounds))

    # Keep track of differences in compounds.
    compoundsOnlyInFirst = list()
    compoundsOnlyInSecond = list()
    diffNames = set()
    diffFormulas = set()
    diffCharges = set()
    different = set()

    # Iterate through all of the compounds in the first object.
    for cpdId in firstCompounds:
        if cpdId not in secondCompounds: # Match on ID
            compoundsOnlyInFirst.append(cpdId)
            continue
        
        # Find compounds with different names.
        if firstCompounds[cpdId]['name'] != secondCompounds[cpdId]['name']:
            diffNames.add(cpdId)
            different.add(cpdId)
            
        # Find compounds with different formulas.
        if firstCompounds[cpdId]['formula'] != secondCompounds[cpdId]['formula']:
            diffFormulas.add(cpdId)
            different.add(cpdId)
        
        # Find compounds with different charges.
        if firstCompounds[cpdId]['charge'] != secondCompounds[cpdId]['charge']:
            diffCharges.add(cpdId)
            different.add(cpdId)
            
    # Find compounds that are only in the second object.
    for cpdId in secondCompounds:
        if cpdId not in firstCompounds:
            compoundsOnlyInSecond.append(cpdId)

    # Print summary data.
    print
    print 'Compound differences:'
    print '  There are %d compound IDs only in first object' %(len(compoundsOnlyInFirst))
    print '  There are %d compound IDs only in second object' %(len(compoundsOnlyInSecond))
    print '  There are %d compounds with same ID and different names' %(len(diffNames))
    print '  There are %d compounds with same ID and different formulas' %(len(diffFormulas))
    print '  There are %d compounds with same ID and different charges' %(len(diffCharges))
    print '  There are %d compounds with same ID and at least one difference' %(len(different))
    diffFormulaCharge = diffFormulas & diffCharges
    print '  There are %d compounds with the same ID and different formulas and charges' %(len(diffFormulaCharge))
    allDiffs = diffNames & diffFormulas & diffCharges
    print '  There are %d compounds with the same ID and all differences' %(len(allDiffs))

    # Print details if requested. 
    if args.showNames:
        if len(diffNames) > 0:
            print 'Compounds with different names:'
            for cpdId in diffNames:
                print '1: line %05d: %s' %(firstCompounds[cpdId]['lineno'], firstCompounds[cpdId])
                print '2: line %05d: %s' %(secondCompounds[cpdId]['lineno'], secondCompounds[cpdId])
                print

    if args.showFormulas:
        if len(diffFormulas) > 0:
            print 'Compounds with different formulas:'
            for cpdId in diffFormulas:
                print '1: line %05d: %s' %(firstCompounds[cpdId]['lineno'], firstCompounds[cpdId])
                print '2: line %05d: %s' %(secondCompounds[cpdId]['lineno'], secondCompounds[cpdId])
                print

    if args.showCharges:
        if len(diffCharges) > 0:
            print 'Compounds with different charges:'
            for cpdId in diffCharges:
                print '1: line %05d: %s' %(firstCompounds[cpdId]['lineno'], firstCompounds[cpdId])
                print '2: line %05d: %s' %(secondCompounds[cpdId]['lineno'], secondCompounds[cpdId])
                print

    exit(0)
