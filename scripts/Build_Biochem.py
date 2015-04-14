#! /usr/bin/python

import argparse
from BiochemHelper import BiochemHelper
from biokbase.workspace.client import Workspace

desc1 = '''
NAME
      Build_Biochem -- build a Biochemistry object

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
'''

desc3 = '''
EXAMPLES
      Build a Biochemistry object:
      > Build_Biochem.py compounds.tsv reactions.tsv
      
SEE ALSO
      Print_Compounds.pl
      Print_Reactions.pl

AUTHORS
      Mike Mundy 
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Build_Biochem', epilog=desc3)
    parser.add_argument('id', help='ID of object', action='store')
    parser.add_argument('cpdfile', help='path to compounds file', action='store')
    parser.add_argument('rxnfile', help='path to reactions file', action='store')
    parser.add_argument('--name', help='name of object', action='store', dest='name', default=None)
    parser.add_argument('--desc', help='description of object', action='store', dest='description', default=None)
    parser.add_argument('-w', '--workspace', help='ID of workspace containing Media object', action='store', dest='workspace', default=None)
    parser.add_argument('--wsurl', help='URL of workspace server', action='store', dest='wsurl', default='https://kbase.us/services/ws/')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # The following fields are required in a Biochemistry object.
    biochem = dict()
    biochem['id'] = args.id
    biochem['compartments'] = list()
    biochem['compounds'] = list()
    biochem['reactions'] = list()
    biochem['reactionSets'] = list()
    biochem['compoundSets'] = list()
    biochem['cues'] = list()
    biochem['compound_aliases'] = dict()
    biochem['reaction_aliases'] = dict()

    # The following fields are optional in a Biochemistry object.
    if args.name is not None:
        biochem['name'] = args.name
    if args.description is not None:
        biochem['description'] = args.description

    # Add the compounds from the compounds file.  Required fields: id, name,
    # name, abbreviation, formula, defaultCharge, isCofactor
    helper = BiochemHelper()
    compounds = helper.readCompoundsFile(args.cpdfile, includeLinenum=False)

    for index in range(len(compounds)):
        biochem['compounds'].append(compounds[index])

    # Add the reactions from the reactions file.  Required fields: id, name,
    # abbreviation, direction, thermoReversibility, status, defaultProtons, reagents.
    reactions = helper.readReactionsFile(args.rxnfile, includeLinenum=False)
    
    for index in range(len(reactions)):
        rxn = reactions[index]
        reactants, products = helper.parseEquation(rxn['equation'])
        if reactants is None and products is None: # @todo Need to confirm this
            continue
        rxn['reagents'] = list()
        for rindex in range(len(reactants)):
            cpd = helper.parseCompoundIdStoich(reactants[rindex])
            reagent = dict()
            reagent['compound_ref'] = '~/compounds/id/'+cpd['compound']
            reagent['compartment_ref'] = '~/compartments/id/'+cpd['compartmentId']
            reagent['coefficient'] = cpd['stoich']*-1.0
            reagent['isCofactor'] = 0 # @todo Is this set separately from value in compound?
            rxn['reagents'].append(reagent)
        for pindex in range(len(products)):
            cpd = helper.parseCompoundIdStoich(products[pindex])
            reagent = dict()
            reagent['compound_ref'] = '~/compounds/id/'+cpd['compound']
            reagent['compartment_ref'] = '~/compartments/id/'+cpd['compartmentId']
            reagent['coefficient'] = cpd['stoich']
            reagent['isCofactor'] = 0 # @todo Is this set separately from value in compound?
            rxn['reagents'].append(reagent)
        del rxn['equation'] # Remove after converting to reagent format
        biochem['reactions'].append(rxn)

    # Save the Biochemistry object to the specified workspace.
    wsClient = Workspace(args.wsurl)
    objectSaveData = dict()
    objectSaveData['type'] = 'KBaseBiochem.Biochemistry-4.0'
    objectSaveData['name'] = args.id
    objectSaveData['data'] = biochem
#    objectSaveData['meta'] = objectMetaData
#    objectSaveData['provenance'] = [ objectProvData ]
    objectInfo = wsClient.save_objects( { 'workspace': args.workspace, 'objects': [ objectSaveData ] } )

    exit(0)