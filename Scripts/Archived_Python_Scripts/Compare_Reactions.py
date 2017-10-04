#! /usr/bin/python

import argparse
from BiochemHelper import BiochemHelper

desc1 = '''
NAME
      Compare_Reactions -- compare two reactions files

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Compare two reactions files and find differences in reactions that have
      have the same ID.  The rxnfile1 and rxnfile2 arguments specify the paths
      to reactions files which describe the reactions in a biochemistry.  The
      first line of the file is a header with the field names.  There is one
      reaction per line with fields separated by tabs.

      When the --show-names optional argument is specified, details on reactions
      with different names are displayed.  When the --show-status optional
      argument is specified, details on reactions with different status are
      displayed.  When the --show-stoich optional argument is specified, details
      on reactions with different stoichiometry are displayed.  When the
      --show-details optional argument is specified, all differences are
      displayed.
'''

desc3 = '''
EXAMPLES
      Show summary data about differences between two reactions files:
      > Compare_Reactions.py reactions1.tsv reactions2.tsv
      
SEE ALSO
      Compare_Compounds

AUTHORS
      Mike Mundy 
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Compare_Reactions', epilog=desc3)
    parser.add_argument('rxnfile1', help='path to first reactions file', action='store')
    parser.add_argument('rxnfile2', help='path to second reactions file', action='store')
    parser.add_argument('--show-details', help='show details on all differences', action='store_true', dest='showDetails', default=False)
    parser.add_argument('--show-names', help='show details on reactions with different names', action='store_true', dest='showNames', default=False)
    parser.add_argument('--show-status', help='show details on reactions with different status', action='store_true', dest='showStatus', default=False)
    parser.add_argument('--show-stoich', help='show details on reactions with different stoichiometry', action='store_true', dest='showStoich', default=False)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # The --show-details option turns on all of the detail options.
    if args.showDetails:
        args.showNames = True
        args.showStatus = True
        args.showStoich = True

    # Read the reactions from the first file.
    helper = BiochemHelper()
    firstReactions = helper.buildDictFromListOfObjects(helper.readReactionsFile(args.rxnfile1))
    print('First reaction file: %s' %(args.rxnfile1))
    print('  Number of reactions: %d' %(len(firstReactions)))

    # Read the compounds from the second file.
    secondReactions = helper.buildDictFromListOfObjects(helper.readReactionsFile(args.rxnfile2))
    print('Second reaction file: %s' %(args.rxnfile2))
    print('  Number of reactions: %d' %(len(secondReactions)))

    # Track differences in reactions.
    reactionsOnlyInFirst = list()
    reactionsOnlyInSecond = list()
    diffNames = set()
    diffStatus = set()
    diffStoich = set()
    different = set()

    # Iterate through all of the reactions in the first object.
    for rxnId in firstReactions:
        if rxnId not in secondReactions: # Match on ID
            reactionsOnlyInFirst.append(rxnId)
            continue
        
        # Find reactions with different names.
        if firstReactions[rxnId]['name'] != secondReactions[rxnId]['name']:
            diffNames.add(rxnId)
            different.add(rxnId)

        # Find reactions with different status.
        if firstReactions[rxnId]['status'] != secondReactions[rxnId]['status']:
            diffStatus.add(rxnId)
            different.add(rxnId)

        # Find reactions with different stoichiometry.  Convert the dictionaries to a tuple
        # so we can use set operations to compare.
        firstReactants = set()
        firstProducts = set()
        reactants, products = helper.parseEquation(firstReactions[rxnId]['equation'])
        for index in range(len(reactants)):
            cpd = helper.parseCompoundIdStoich(reactants[index])
            firstReactants.add( ( cpd['stoich'], cpd['compound'] ) )
        for index in range(len(products)):
            cpd = helper.parseCompoundIdStoich(products[index])
            firstProducts.add( ( cpd['stoich'], cpd['compound'] ) )
        secondReactants = set()
        secondProducts = set()
        reactants, products = helper.parseEquation(secondReactions[rxnId]['equation'])
        for index in range(len(reactants)):
            cpd = helper.parseCompoundIdStoich(reactants[index])
            secondReactants.add( ( cpd['stoich'], cpd['compound'] ) )
        for index in range(len(products)):
            cpd = helper.parseCompoundIdStoich(products[index])
            secondProducts.add( ( cpd['stoich'], cpd['compound'] ) )
        diffReactants = firstReactants ^ secondReactants # Elements in either set but not both
        diffProducts = firstProducts ^ secondProducts
        if len(diffReactants) > 0 or len(diffProducts) > 0:
            diffStoich.add(rxnId)
            different.add(rxnId)

    # Find reactions that are only in the second object.
    for rxnId in secondReactions:
        if rxnId not in firstReactions:
            reactionsOnlyInSecond.append(rxnId)

    # Print summary data.
    print()
    print('Reaction differences:')
    print('  There are %d reaction IDs only in first object' %(len(reactionsOnlyInFirst)))
    print('  There are %d reaction IDs only in second object' %(len(reactionsOnlyInSecond)))
    print('  There are %d reactions with same ID and different names' %(len(diffNames)))
    print('  There are %d reactions with same ID and different status' %(len(diffStatus)))
    print('  There are %d reactions with same ID and different stoichiometry' %(len(diffStoich)))
    print('  There are %d reactions with same ID and at least one difference' %(len(different)))
    allDiffs = diffNames & diffStatus & diffStoich
    print('  There are %d reactions with the same ID and all differences' %(len(allDiffs)))
    print()
    
    # Print details if requested. 
    if args.showNames:
        if len(diffNames) > 0:
            print('Reactions with different names:')
            for rxnId in diffNames:
                print('1: line %05d: %s' %(firstReactions[rxnId]['linenum'], firstReactions[rxnId]))
                print('2: line %05d: %s' %(secondReactions[rxnId]['linenum'], secondReactions[rxnId]))
                print()

    if args.showStatus:
        if len(diffStatus) > 0:
            print('Reactions with different status:')
            for rxnId in diffStatus:
                print('1: line %05d: %s' %(firstReactions[rxnId]['linenum'], firstReactions[rxnId]))
                print('2: line %05d: %s' %(secondReactions[rxnId]['linenum'], secondReactions[rxnId]))
                print()

    if args.showStoich:
        if len(diffStoich) > 0:
            print('Reactions with different stoichiometry:')
            for rxnId in diffStoich:
                print('1: line %05d: %s' %(firstReactions[rxnId]['linenum'], firstReactions[rxnId]))
                print('2: line %05d: %s' %(secondReactions[rxnId]['linenum'], secondReactions[rxnId]))
                print()

    exit(0)
