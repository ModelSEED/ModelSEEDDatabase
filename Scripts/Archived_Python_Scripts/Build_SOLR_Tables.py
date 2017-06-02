#! /usr/bin/env python

import argparse
import os
from BiochemHelper import BiochemHelper

desc1 = '''
NAME
      Build_SOLR_Tables -- build tables for importing to SOLR

SYNOPSIS
      Build the tables that are used to import the biochemistry data into SOLR.
'''

desc2 = '''
DESCRIPTION
'''

desc3 = '''
EXAMPLES
      Build SOLR tables:
      > Build_SOLR_Tables.py --outputdir ../SOLRTables
      
SEE ALSO
      Build_Biochem.py

AUTHORS
      Mike Mundy 
'''

if __name__ == '__main__':
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Build_SOLR_Tables', epilog=desc3)
    parser.add_argument('--outputdir', help='path to output directory for storing SOLR tables', action='store', default='../SOLRDump')
    parser.add_argument('--compoundfile', help='path to source master compounds file', action='store', default='../Biochemistry/compounds.master.tsv')
    parser.add_argument('--reactionfile', help='path to source master reactions file', action='store', default='../Biochemistry/reactions.master.tsv')
    parser.add_argument('--aliasdir', help='path to directory with source aliases files', action='store', default='../Aliases')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Create a helper for working with source files.
    helper = BiochemHelper()

    # Get the aliases from all of the aliases files.
    compoundAliases, reactionAliases = helper.readAliasFiles(args.aliasdir)

    # Get the compounds from the master compounds file.
    compounds = helper.readCompoundsFile(args.compoundfile, includeLinenum=False, noFormat=True)
    
    # Build the compounds file.
    compoundHeader = [ 'id', 'abbreviation', 'name', 'formula', 'mass', 'source',
                       'structure', 'charge', 'is_core', 'is_obsolete', 'linked_compound',
                       'is_cofactor', 'deltag', 'deltagerr', 'pka', 'pkb',
                       'abstract_compound', 'comprised_of', 'aliases' ]

    compoundFile = os.path.join(args.outputdir, 'Compounds.tsv')
    with open(compoundFile, 'w') as handle:
        handle.write('%s\n' %('\t'.join(compoundHeader)))
        for index in range(len(compounds)):
            compound = compounds[index]
            if compound['id'] in compoundAliases:
                aliasList = list()
                for source in compoundAliases[compound['id']]:
                    for aindex in range(len(compoundAliases[compound['id']][source])):
                        aliasList.append('"%s:%s"' %(source, compoundAliases[compound['id']][source][aindex]))
                compound['aliases'] = ';'.join(aliasList)
                line = [ compound[h] for h in compoundHeader ]
                handle.write('%s\n' %('\t'.join(line)))
    print 'Stored compounds in '+compoundFile

    # Get the reactions from the master reactions file.
    reactions = helper.readReactionsFile(args.reactionfile, includeLinenum=False, noFormat=True)

    # Build the reactions file.
    reactionHeader = [ 'id', 'abbreviation', 'name', 'code', 'stoichiometry', 'is_transport', 
                       'equation', 'definition', 'reversibility', 'direction', 'abstract_reaction',
                       'pathways', 'aliases', 'ec_numbers', 'deltag', 'deltagerr', 'compound_ids',
                       'status' ]
    
    reactionFile = os.path.join(args.outputdir, 'Reactions.tsv')
    with open(reactionFile, 'w') as handle:
        handle.write('%s\n' %('\t'.join(reactionHeader)))
        for index in range(len(reactions)):
            reaction = reactions[index]
            if reaction['id'] in reactionAliases:
                aliasList = list()
                for source in reactionAliases[reaction['id']]:
                    for aindex in range(len(reactionAliases[reaction['id']][source])):
                        aliasList.append('"%s:%s"' %(source, reactionAliases[reaction['id']][source][aindex]))
                reaction['aliases'] = ';'.join(aliasList)
                line = [ reaction[h] for h in reactionHeader ]
                handle.write('%s\n' %('\t'.join(line)))
    print 'Stored reactions in '+reactionFile
    
    exit(0)
