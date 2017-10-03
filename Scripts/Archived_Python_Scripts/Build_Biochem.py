#! /usr/bin/env python

import argparse
import json
from BiochemHelper import BiochemHelper
from biop3.Workspace.WorkspaceClient import Workspace

desc1 = '''
NAME
      Build_Biochem -- build a Biochemistry typed object

SYNOPSIS
      Build a Biochemistry typed object directly from source files.
'''

desc2 = '''
DESCRIPTION
'''

desc3 = '''
EXAMPLES
      Build a Biochemistry object:
      > Build_Biochem.py master-2015a /mmundy/public/modelsupport/biochemistry/master-2015a
      
SEE ALSO
      Build_SOLR_Tables.py

AUTHORS
      Mike Mundy 
'''

if __name__ == '__main__':
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Build_Biochem', epilog=desc3)
    parser.add_argument('id', help='ID of Biochemistry object', action='store')
    parser.add_argument('ref', help='reference to workspace location to store Biochemistry object')
#    parser.add_argument('compartmentfile', help='path to compartments file', action='store')
    parser.add_argument('--compoundfile', help='path to compounds file', action='store', default='../Biochemistry/compounds.master.tsv')
    parser.add_argument('--reactionfile', help='path to reactions file', action='store', default='../Biochemistry/reactions.master.tsv')
    parser.add_argument('--aliasdir', help='path to directory with source aliases files', action='store', default='../Aliases')
    parser.add_argument('--name', help='name of object', action='store', dest='name', default=None)
    parser.add_argument('--desc', help='description of object', action='store', dest='description', default=None)
    parser.add_argument('--wsurl', help='URL of workspace server', action='store', dest='wsurl', default='http://p3.theseed.org/services/Workspace')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Create the Biochemistry typed object.
    biochem = dict()
    biochem['id'] = args.id

    # The following fields are required in a Biochemistry typed object but they are currently unused.
    biochem['reactionSets'] = list()
    biochem['compoundSets'] = list()
    biochem['cues'] = list()

    # The following fields are optional in a Biochemistry typed object.
    if args.name is not None:
        biochem['name'] = args.name
    if args.description is not None:
        biochem['description'] = args.description

    # Create a helper object.
    helper = BiochemHelper()

    # Add the compounds from the compounds file.
    print(('Adding compounds from %s ...' %(args.compoundfile)))
    biochem['compounds'] = helper.readCompoundsFile(args.compoundfile, includeLinenum=False)
    compounds = helper.buildIndexDictFromListOfObjects(biochem['compounds'])

    # Start with an empty dictionary of compartments.  With compartment-free reactions, just add
    # place holder compartments as they are found processing the reactions.
    compartments = dict()

    # Add the reactions from the reactions file.
    print(('Adding reactions from %s ...' %(args.reactionfile)))
    biochem['reactions'] = list()
    reactions = helper.readReactionsFile(args.reactionfile, includeLinenum=False)
    
    for index in range(len(reactions)):
        rxn = reactions[index]
        reactants, products = helper.parseEquation(rxn['equation'])
        if reactants is None and products is None: # @todo Need to confirm this
            continue
        rxn['reagents'] = list()
        for rindex in range(len(reactants)):
            cpd = helper.parseCompoundIdStoich(reactants[rindex])
            reagent = dict()
            # Validate the compound is valid.
            if cpd['compound'] in compounds:
                reagent['compound_ref'] = '~/compounds/id/'+cpd['compound']
            else:
                print(('WARNING: Compound %s is not defined in the list of compounds' %(cpd['compound'])))
            # Add compartment the first time it is found.
            if cpd['compartmentId'] not in compartments:
                compartments[cpd['compartmentId']] = { 'id': cpd['compartmentId'], 'name': 'Compartment'+cpd['compartmentId'], 'hierarchy': 3}
            reagent['compartment_ref'] = '~/compartments/id/'+cpd['compartmentId']
            reagent['coefficient'] = cpd['stoich']*-1.0
            reagent['isCofactor'] = 0 # @todo Is this set separately from value in compound?
            rxn['reagents'].append(reagent)
        for pindex in range(len(products)):
            cpd = helper.parseCompoundIdStoich(products[pindex])
            reagent = dict()
            # Validate the compound is valid.
            if cpd['compound'] in compounds:
                reagent['compound_ref'] = '~/compounds/id/'+cpd['compound']
            else:
                print(('WARNING: Compound %s is not defined in the list of compounds' %(cpd['compound'])))
            # Add compartment the first time it is found.
            if cpd['compartmentId'] not in compartments:
                compartments[cpd['compartmentId']] = { 'id': cpd['compartmentId'], 'name': 'Compartment'+cpd['compartmentId'], 'hierarchy': 3}
            reagent['compartment_ref'] = '~/compartments/id/'+cpd['compartmentId']
            reagent['coefficient'] = cpd['stoich']
            reagent['isCofactor'] = 0 # @todo Is this set separately from value in compound?
            rxn['reagents'].append(reagent)
        del rxn['equation'] # Remove after converting to reagent format
        biochem['reactions'].append(rxn)

    # Create the compartment list from the dictionary assembled above.
    biochem['compartments'] = list()
    for id in compartments:
        biochem['compartments'].append(compartments[id])

    # Add the aliases from all of the aliases files.
    print(('Reading aliases from %s ...' %(args.aliasdir)))
    compoundAliases, reactionAliases = helper.readAliasFiles(args.aliasdir)
    biochem['compound_aliases'] = compoundAliases
    biochem['reaction_aliases'] = reactionAliases

    # Save the Biochemistry typed object to the specified workspace path. An existing typed object
    # is overwritten with the updated data.
    print(('Saving typed object to %s ...' %(args.ref)))
    wsClient = Workspace(args.wsurl)
    output = wsClient.create( { 'objects': [ [ args.ref, 'biochemistry', {}, biochem ] ], 'overwrite': 1 });
    
#    json.dump(biochem, open('master.json', 'w'), indent=4)
    exit(0)