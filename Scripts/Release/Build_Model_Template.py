#! /usr/bin/env python

import argparse
import os
import json
from Scripts.TemplateHelper import TemplateHelper

desc1 = '''
NAME
      Build_Model_Template -- build a Model Template object from source files

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
'''

desc3 = '''
EXAMPLES
      Build a Model Template object:
      > Build_Model_Template.py GramNegative ../Templates/GramNegative /mmundy/public/modelsupport/templates/GramNegative.modeltemplate
      
SEE ALSO
      Build_Biochem.py

AUTHORS
      Mike Mundy 
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Build_Model_Template', epilog=desc3)
    parser.add_argument('id', help='ID of Model Template object', action='store')
    parser.add_argument('templatedir', help='path to directory containing source files', action='store')
    parser.add_argument('--biochemref', help='reference to Biochemistry object in workspace', action='store', default='default/default.biochem')
    parser.add_argument('--compoundfile', help='path to master compounds file', action='store', default='../../Biochemistry/compounds.json')
    parser.add_argument('--reactionfile', help='path to master reactions file', action='store', default='../../Biochemistry/reactions.json')
    parser.add_argument('--complexfile', help='path to master complexes file', action='store', default='../../Annotations/Complexes.tsv')
    parser.add_argument('--rolefile', help='path to master roles file', action='store', default='../../Annotations/Roles.tsv')
    parser.add_argument('--name', help='name of object', action='store', default=None)
    parser.add_argument('--type', help='type of model', action='store', default='GenomeScale')
    parser.add_argument('--domain', help='domain of organisms', action='store', default='Bacteria')
    parser.add_argument('--show-stats', help='show statistics about Model Template', action='store_true', dest='showStats', default=False)
    parser.add_argument('--wsurl', help='URL of workspace server', action='store', dest='wsurl', default='http://p3.theseed.org/services/Workspace')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Create a helper object.
    helper = TemplateHelper(args.compoundfile, args.reactionfile)

    # The following fields are required in a Model Template object.
    template = dict()
    template['id'] = args.id
    if args.name is not None:
        template['name'] = args.name
    else:
        template['name'] = args.id
    template['type'] = args.type
    template['domain'] = args.domain
    template['biochemistry_ref'] = args.biochemref
    template['pathways'] = list() # Always an empty for now

    # Order is important so references can be made between sections of the Model Template.
    
    # Add the template compartments.
    compartmentsFile = os.path.join(args.templatedir, 'Compartments.tsv')
    helper.readCompartmentsFile(compartmentsFile, includeLinenum=False)
    template['compartments'] = [ helper.compartments[key] for key in helper.compartments ]

    # Add the template biomasses.
    biomassFile = os.path.join(args.templatedir, 'Biomasses.tsv')
    biomassCompoundsFile = os.path.join(args.templatedir, 'BiomassCompounds.tsv')
    helper.readBiomassesFile(biomassFile, biomassCompoundsFile, includeLinenum=False)
    template['biomasses'] = [ helper.biomasses[key] for key in helper.biomasses ]

    # Add the template roles.
    helper.readRolesFile(args.rolefile, includeLinenum=False)
    template['roles'] = [ helper.roles[key] for key in helper.roles ]

    # Add the template complexes.
    helper.readComplexesFile(args.complexfile, includeLinenum=False)
    template['complexes'] = [ helper.complexes[key] for key in helper.complexes ]

    # Add the template reactions.
    reactionsFile = os.path.join(args.templatedir, 'Reactions.tsv')
    helper.readReactionsFile(reactionsFile, includeLinenum=False)
    template['reactions'] = [ helper.reactions[key] for key in helper.reactions ]

    # Add the template compounds (constructed from reagents in reactions).
    template['compounds'] = [ helper.compounds[key] for key in helper.compounds ]

    # Add the template comp compounds (constructed from reagents in reactions).
    template['compcompounds'] = [ helper.compCompounds[key] for key in helper.compCompounds ]

    # Save a local copy for easy reference.
    filename = os.path.join(args.templatedir, args.id+'.json')
    print(filename)
    json.dump(template, open(filename, 'w'), indent=4)
    
    # Print statistics about the Model Template if requested.
    if args.showStats:
        print(('Number of compartments: '+str(len(template['compartments']))))
        print(('Number of biomasses: '+str(len(template['biomasses']))))
        print(('Number of roles: '+str(len(template['roles']))))
        print(('Number of complexes: '+str(len(template['complexes']))))
        print(('Number of reactions: '+str(len(template['reactions']))))
        print(('  Number of conditional reactions: '+str(helper.numConditional)))
        print(('  Number of gapfilling reactions: '+str(helper.numGapfilling)))
        print(('  Number of spontaneous reactions: '+str(helper.numSpontaneous)))
        print(('  Number of universal reactions: '+str(helper.numUniversal)))
        print(('Number of compounds: '+str(len(template['compounds']))))
        print(('Number of compcompounds: '+str(len(template['compcompounds']))))
        
    exit(0)
