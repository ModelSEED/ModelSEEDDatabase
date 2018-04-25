#! /usr/bin/env python

import argparse
import re
import sys
from collections import Counter
from csv import DictReader

from repostat.stash import StatStash
from .error_reporting import find_new_errors, report_errors
from ..Biochem_Helper import BiochemHelper

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


def get_atom_count(formulaDict, complist):
    atom_counts = Counter()
    for term in complist:
        stoich = float(term.split()[0].strip("()"))
        id = term.split()[1].split('[')[0]
        if formulaDict[id] == 'null':
            continue
        for pair in re.findall('([A-Z][a-z]?)(\d*)', formulaDict[id]):
            if not pair[1]:
                atom_counts[pair[0]] += stoich
            else:
                atom_counts[pair[0]] += int(pair[1]) * stoich
    return atom_counts


if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Validate_Reactions', epilog=desc3)
    parser.add_argument('rxnfile', help='path to reactions file', action='store')
    parser.add_argument('-c', help='Path to compound file', action='store', dest='compfile', default=False)
    parser.add_argument('--show-details', help='show details on all problems', action='store_true', dest='showDetails', default=False)
    parser.add_argument('--show-dup-ids', help='show details on duplicate IDs', action='store_true', dest='showDupIds', default=False)
    parser.add_argument('--show-bad-ids', help='show details on bad IDs', action='store_true', dest='showBadIds', default=False)
    parser.add_argument('--show-dup-names', help='show details on duplicate names', action='store_true', dest='showDupNames', default=False)
    parser.add_argument('--show-bad-names', help='show details on bad names', action='store_true', dest='showBadNames', default=False)
    parser.add_argument('--show-dup-abbrs', help='show details on duplicate abbreviations', action='store_true', dest='showDupAbbrs', default=False)
    parser.add_argument('--show-bad-abbrs', help='show details on bad abbreviations', action='store_true', dest='showBadAbbrs', default=False)
    parser.add_argument('--show-bad-direction', help='show details on bad directions', action='store_true', dest='showBadDirection', default=False)
    parser.add_argument('--show-bad-reverse', help='show details on bad reversibility', action='store_true', dest='showBadReverse', default=False)
    parser.add_argument('--show-diff-eq', help='show details on different equation and code', action='store_true', dest='showDiffEqCode', default=False)
    parser.add_argument('--show-dup-eq', help='show details on duplicated reaction equations', action='store_true', dest='showDupEquations', default=False)
    parser.add_argument('--show-unbalanced', help='show details on equations taht are not atom balanced', action='store_true', dest='showUnbalanced', default=False)
    parser.add_argument('--show-bad-eq', help='show details on missing reactants or products in equations', action='store_true', dest='showBadEquations', default=False)
    parser.add_argument('--show-status', help='show details on status types', action='store_true', dest='showStatus', default=False)
    parser.add_argument('--show-bad-link', help='show details on bad links', action='store_true', dest='showBadLink', default=False)
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
        args.showBadAbbrs = True
        args.showDupEquations = True
        args.showStatus = True
        args.showUnbalanced = True

    # Read the reactions from the specified file.
    print('Reaction file: %s' % args.rxnfile)
    helper = BiochemHelper()
    reactions = helper.readReactionsFile(args.rxnfile)
    if reactions is None:
        print('Error reading reactions file')
        exit(1)
    print('Number of reactions: %d' % len(reactions))
    
    # Create a dictionary keyed by id for fast lookup of reactions (and compounds)
    reactionDict = helper.buildIndexDictFromListOfObjects(reactions)
    if args.compfile:
        compoundDict = {}
        obs_compounds = set()
        with open(args.compfile, 'r') as infile:
            for line in DictReader(infile, dialect='excel-tab'):
                compoundDict[line['id']] = line['formula']
                if line['is_obsolete']:
                    obs_compounds.add(line['id'])

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
    badDirection = list()
    badReversibility = list()
    unknownReversibility = list()
    diffEquationCode = list()
    noEquation = list()
    duplicateEquation = dict()
    unbalanced = list()
    obsoleteComps = list()
    eqnHashDict = dict()
    noDefinition = list()
    noReactants = list()
    noProducts = list()
    statusTypes = {'CK': 0, 'OK': 0, 'MI': 0, 'CI': 0, 'HB': 0, 'EMPTY': 0,
                   'CPDFORMERROR': 0}
    okStatus = 0
    badTransport = list()
    isTransport = list()
    badObsolete = list()
    isObsolete = list()
    unknownDeltag = list()
    unknownDeltagErr = list()
    zeroDeltag = list()
    zeroDeltagerr = list()
    badLink = list()
    
    for index in range(len(reactions)):
        rxn = reactions[index]

        # Check for invalid is_obsolete flags.
        if 'is_obsolete' in rxn:
            if rxn['is_obsolete'] != 0 and rxn['is_obsolete'] != 1:
                badObsolete.append(index)
            if rxn['is_obsolete'] == 1:
                isObsolete.append(index)
                continue
        
        # Check for duplicate IDs.
        if rxn['id'] in idDict:
            if len(idDict[rxn['id']]) == 1:
                duplicateId += 1
            idDict[rxn['id']].append(index)
        else:
            idDict[rxn['id']] = [index]

        # Check for invalid characters in the ID.
        match = re.search(r'^rxn\d\d\d\d\d$', rxn['id'])
        if match is None:
            badIdChars.append(index)

        # Check for duplicate names.
        if rxn['name'] in nameDict:
            if len(nameDict[rxn['name']]) == 1:
                duplicateName += 1
            nameDict[rxn['name']].append(index)
        else:
            nameDict[rxn['name']] = [index]

        # Check for invalid characters in the name.
        try:
            rxn['name'].encode('ascii')
        except UnicodeEncodeError:
            badNameChars.append(index)
        if rxn['name'] != rxn['name'].strip():
            badNameChars.append(index)

        # Check for duplicate abbreviations.
        if rxn['abbreviation'] in abbrDict:
            if len(abbrDict[rxn['abbreviation']]) == 1:
                duplicateAbbr += 1
            abbrDict[rxn['abbreviation']].append(index)
        else:
            abbrDict[rxn['abbreviation']] = [index]

        # Check for invalid characters in the abbreviation.
        try:
            rxn['abbreviation'].encode('ascii')
        except UnicodeEncodeError:
            badAbbrChars.append(index)
        if rxn['abbreviation'] != rxn['abbreviation'].strip():
            badAbbrChars.append(index)

        # Check for invalid direction.
        if rxn['direction'] != '<' and rxn['direction'] != '>' and rxn['direction'] != '=':
            badDirection.append(index)

        # Check for unknown or invalid reversibility.
        if rxn['reversibility'] == '?':
            unknownReversibility.append(index)
        elif rxn['reversibility'] != '<' and rxn['reversibility'] != '>' and rxn['reversibility'] != '=':
            badReversibility.append(index)

        # Check for different equation and code fields.
        if rxn['equation'] != rxn['code']:
            diffEquationCode.append(index)

        # Check for missing reactants and/or products.
        reactants, products = helper.parseEquation(rxn['equation'])
        if args.compfile:
            reactant_atoms = get_atom_count(compoundDict, reactants)
            product_atoms = get_atom_count(compoundDict, products)
            if reactant_atoms - product_atoms or product_atoms - reactant_atoms:
                unbalanced.append((index, reactant_atoms, product_atoms))
            if [cid for cid in reactants + products if cid in obs_compounds]:
                obsoleteComps.append(index)

        if reactants is None and products is None:
            noEquation.append(index)
        else:
            rxn_hash = hash("%s<>%s" % (sorted(reactants), sorted(products)))
            if rxn_hash in eqnHashDict:
                if rxn_hash not in duplicateEquation:
                    duplicateEquation[rxn_hash] = [eqnHashDict[rxn_hash]]
                duplicateEquation[rxn_hash].append(index)
            else:
                eqnHashDict[rxn_hash] = index

            if len(reactants) == 0:
                noReactants.append(index)
            if len(products) == 0:
                noProducts.append(index)

        # Check reaction status.
        if rxn['status'] == 'OK':
            okStatus += 1
        else:
            fields = rxn['status'].split('|')
            simpleStatus = list()
            for sindex in range(len(fields)):
                if ':' in fields[sindex]:
                    pos = fields[sindex].find(':')
                    type = fields[sindex][:pos]
                else:
                    type = fields[sindex]
                statusTypes[type] += 1
                simpleStatus.append(type)
            if len(simpleStatus) > 1:
                status = '|'.join(simpleStatus)
                if status in statusTypes:
                    statusTypes[status] += 1
                else:
                    statusTypes[status] = 1

        # Check for invalid is_transport flags.
        if rxn['is_transport'] != 0 and rxn['is_transport'] != 1:
            badTransport.append(index)
        if rxn['is_transport'] == 1:
            isTransport.append(index)

        # Check that linked reactions are all valid.
        if 'linked_reaction' in rxn:
            linkedRxns = rxn['linked_reaction'].split(';')
            for rxnid in linkedRxns:
                if rxnid not in reactionDict:
                    badLink.append(index)
        
        # Check for unknown deltaG and deltaGerr values.
        if 'deltag' in rxn:
            if rxn['deltag'] == float(0):
                zeroDeltag.append(index)
        else:
            unknownDeltag.append(index)
        if 'deltagerr' in rxn:
            if rxn['deltagerr'] == float(0):
                zeroDeltagerr.append(index)
        else:
            unknownDeltagErr.append(index)

    # Print summary data.
    print('Number of reactions with duplicate IDs: %d' % (duplicateId))
    print('Number of reactions with bad characters in ID: %d' % (len(badIdChars)))
    print('Number of reactions with duplicate names: %d' % (duplicateName))
    print('Number of reactions with bad characters in name: %d' % (len(badNameChars)))
    print('Number of reactions with duplicate abbreviations: %d' % (duplicateAbbr))
    print('Number of reactions with bad characters in abbreviation: %d' % (len(badAbbrChars)))
    print('Number of reactions with bad direction: %d' % (len(badDirection)))
    print('Number of reactions with bad reversibility: %d' % (len(badReversibility)))
    print('Number of reactions with unknown reversibility: %d' % (len(unknownReversibility)))
    print('Number of reactions with different equation and code: %d' % (len(diffEquationCode)))
    print('Number of reactions with missing equation: %d' % (len(noEquation)))
    print('Number of reactions with duplicate equations: %d' % (len(duplicateEquation)))
    print('Number of reactions with unbalanced equations: %d' % (len(unbalanced)))
    print('Number of reactions with obsolete compounds: %d' % (len(obsoleteComps)))
    print('Number of reactions with no reactants: %d' % (len(noReactants)))
    print('Number of reactions with no products: %d' % (len(noProducts)))
    print('Number of reactions with OK status: %d' % (okStatus))
    print('Number of reactions with OK status after balancing: %d' % (statusTypes['OK']))
    print('Number of reactions with bad status: %d' % (len(reactions)-okStatus-statusTypes['OK']))
    print('Number of reactions with bad is_transport flag: %d' % (len(badTransport)))
    print('Number of transport reactions: %d' % (len(isTransport)))
    print('Number of reactions with bad is_obsolete flag: %d' % (len(badObsolete)))
    print('Number of obsolete reactions: %d' % (len(isObsolete)))
    print('Number of reactions with unknown deltaG value: %d' % (len(unknownDeltag)))
    print('Number of reactions with zero deltaG value: %d' % (len(zeroDeltag)))
    print('Number of reactions with unknown deltaGErr value: %d' % (len(unknownDeltagErr)))
    print('Number of reactions with zero deltaGErr value: %d' % (len(zeroDeltagerr)))
    print('Number of reactions with bad links: %d' % (len(badLink)))
    print('\n')

    # Print details if requested.
    if args.showDupIds:
        for id in idDict:
            if len(idDict[id]) > 1:
                print('Duplicate reaction ID: %s' % id)
                for dup in idDict[id]:
                    print('Line %05d: %s' % (reactions[dup]['linenum'], reactions[dup]))
                print()
    if args.showBadIds:
        if len(badIdChars) > 0:
            print('reactions with bad characters in ID:')
            for index in range(len(badIdChars)):
                print('Line %05d: %s' % (reactions[badIdChars[index]]['linenum'], reactions[badIdChars[index]]))
            print()
    if args.showDupNames:
        for name in nameDict:
            if len(nameDict[name]) > 1:
                print('Duplicate reaction name: %s' % name)
                for dup in nameDict[name]:
                    print('Line %05d: %s' % (reactions[dup]['linenum'], reactions[dup]['id']))
                print()
    if args.showBadNames:
        if len(badNameChars) > 0:
            print('reactions with bad characters in name:')
            for index in range(len(badNameChars)):
                print('Line %05d: %s' % (reactions[badNameChars[index]]['linenum'], reactions[badNameChars[index]]))
            print()
    if args.showDupAbbrs:
        for abbr in abbrDict:
            if len(abbrDict[abbr]) > 1:
                print('Duplicate reaction abbreviation: %s' % abbr)
                for dup in abbrDict[abbr]:
                    print('Line %05d: %s' % (reactions[dup]['linenum'], reactions[dup]))
                print()
    if args.showBadAbbrs:
        if len(badAbbrChars) > 0:
            print('Reactions with bad characters in abbreviation:')
            for index in range(len(badAbbrChars)):
                print('Line %05d: %s' % (reactions[badAbbrChars[index]]['linenum'], reactions[badAbbrChars[index]]))
            print()
    if args.showBadDirection:
        if len(badDirection) > 0:
            print('Reactions with bad value in direction:')
            for index in range(len(badDirection)):
                print('Line %05d: %s' % (reactions[badDirection[index]]['linenum'], reactions[badDirection[index]]))
            print()
    if args.showBadReverse:
        if len(badReversibility) > 0:
            print('Reactions with bad value in reversibility:')
            for index in range(len(badReversibility)):
                print('Line %05d: %s' % (reactions[badReversibility[index]]['linenum'], reactions[badReversibility[index]]))
            print()
        if len(unknownReversibility) > 0:
            print('Reactions with unknown reversibility:')
            for index in range(len(unknownReversibility)):
                print('Line %05d: %s' % (reactions[unknownReversibility[index]]['linenum'], reactions[unknownReversibility[index]]))
            print()
    if args.showDiffEqCode:
        if len(diffEquationCode) > 0:
            print('Reactions with different equation and code fields:')
            for index in range(len(diffEquationCode)):
                print('Line %05d: %s' % (reactions[diffEquationCode[index]]['linenum'], reactions[diffEquationCode[index]]))
            print()
    if args.showBadEquations:
        if len(noReactants) > 0:
            print('Reactions with no reactants:')
            for index in range(len(noReactants)):
                print('Line %05d: %s' % (reactions[noReactants[index]]['linenum'], reactions[noReactants[index]]))
            print()
        if len(noProducts) > 0:
            print('Reactions with no products:')
            for index in range(len(noProducts)):
                print('Line %05d: %s' % (reactions[noProducts[index]]['linenum'], reactions[noProducts[index]]))
            print()
        if len(noEquation) > 0:
            print('Reactions with no equation:')
            for index in range(len(noEquation)):
                print('Line %05d: %s' % (reactions[noEquation[index]]['linenum'], reactions[noEquation[index]]))
            print()
    if args.showDupEquations:
        for rxn_list in duplicateEquation.values():
            print('\nReactions with duplicate equations:')
            for index in rxn_list:
                print('Line %05d: %s' % (reactions[index]['linenum'], reactions[index]['id']))
    if args.showUnbalanced:
        if len(unbalanced) > 0:
            print('Unbalenced Reactions:')
            for tup in unbalanced:
                print('Line %05d: %s' % (reactions[tup[0]]['linenum'], reactions[tup[0]]['id']))
                print("Reactant atoms:%s\n Product atoms:%s\n" %
                      (sorted(tup[1].items()), sorted(tup[2].items())))
    if args.showStatus:
        print('Reactions with status that is not OK:')
        for type in statusTypes:
            print('Reaction status %s: %d' % (type, statusTypes[type]))
    if args.showBadLink:
        if len(badLink) > 0:
            print('Reactions with bad links:')
            for index in range(len(badLink)):
                print('Line %05d: %s' % (reactions[badLink[index]]['linenum'], reactions[badLink[index]]))

    error_fields = ['duplicateId', 'duplicateName', 'duplicateEquation',
                    'unbalanced', 'obsoleteComps', 'badIdChars', 'badLink',
                    'badNameChars', 'badAbbrChars', 'badDirection',
                    'badReversibility', 'badObsolete', 'badTransport']
    errors = dict([(x, eval(x)) if isinstance(eval(x), int)
                   else (x, len(eval(x))) for x in error_fields])
    new_errors = find_new_errors('reactions', errors)
    report_errors('reactions', errors)
    stash = StatStash('redis://redis-16221.c12.us-east-1-4.ec2.cloud.redislabs.com:16221')
    stash.report_stats('reactions', errors)


    if new_errors:
        print("ERRORS: " + ", ".join(new_errors), file=sys.stderr)
        exit(1)
