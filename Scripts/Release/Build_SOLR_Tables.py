#! /usr/bin/env python

import argparse
import os
import csv
import json
from collections import defaultdict
from Scripts.Biochem_Helper import BiochemHelper

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


def alias_dict(path):
    alias_mappings = defaultdict(list)
    with open(path) as infile:
        r = csv.DictReader(infile, dialect='excel-tab')
        for line in r:
            for seed_id in line['MS ID'].split('|'):
                alias_mappings[seed_id].append('"%s:%s"' % (
                    line['Source'].strip(), line['External ID']))
    return alias_mappings


def make_tsv(data, tsv_file, tsv_header):
    with open(tsv_file, 'w') as handle:
        handle.write('%s\n' % ('\t'.join(tsv_header)))
        for item in data:
            line = [str(item[h]) for h in tsv_header]
            handle.write('%s\n' % ('\t'.join(line)))


if __name__ == '__main__':
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Build_SOLR_Tables', epilog=desc3)
    parser.add_argument('--outputdir', help='path to output directory for storing SOLR tables', action='store', default='../../Objects')
    parser.add_argument('--compoundfile', help='path to source master compounds file', action='store', default='../../Biochemistry/compounds.tsv')
    parser.add_argument('--reactionfile', help='path to source master reactions file', action='store', default='../../Biochemistry/reactions.tsv')
    parser.add_argument('--aliasdir', help='path to directory with source aliases files', action='store', default='../../Biochemistry/Aliases')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Create a helper for working with source files.
    helper = BiochemHelper()

    # Get the aliases from all of the aliases files.
    cpd_aliases = alias_dict(args.aliasdir+'/Compounds_Aliases.tsv')
    rxn_aliases  = alias_dict(args.aliasdir+'/Reactions_Aliases.tsv')

    # Get additional cpd_names
    cpd_names = alias_dict(args.aliasdir+'/Names_Compounds_Aliases.tsv')

    # Get the compounds from the master compounds file.
    compounds = helper.readCompoundsFile(args.compoundfile, includeLinenum=False, noFormat=True)
    
    # Build the compounds file.
    compound_header = ['id', 'abbreviation', 'name', 'formula', 'mass',
                       'source', 'structure', 'charge', 'is_core',
                       'is_obsolete', 'is_cofactor', 'deltag', 'deltagerr',
                       'pka', 'pkb', 'aliases']
    
    # Add aliases
    for item in compounds:
        item['aliases'] = ";".join(cpd_aliases[item['id']])

    # Add names
    for item in compounds:
        names = ";".join(filter(lambda x: ("searchname" not in x), cpd_names[item['id']]))
        item['aliases']+=";"+names

    # Write TSV
    compound_file = os.path.join(args.outputdir, 'Compounds.tsv')
    make_tsv(compounds, compound_file, compound_header)

    # Write JSON
    compounds_json = open(os.path.join(args.outputdir, 'Compounds.json'), 'w')
    json.dump(compounds, compounds_json, sort_keys=True, indent=4)

    print('Stored compounds in ' + compound_file + ' and Compounds.json')

    # Get the reactions from the master reactions file.
    reactions = helper.readReactionsFile(args.reactionfile, includeLinenum=False, noFormat=True)

    # Build the reactions file.
    reaction_header = ['id', 'name', 'code', 'stoichiometry', 'is_transport',
                       'equation', 'definition', 'reversibility', 'direction',
                       'aliases', 'deltag', 'deltagerr', 'compound_ids']

    # Add aliases
    for item in reactions:
        item['aliases'] = ";".join(rxn_aliases[item['id']])

    # Write TSV    
    reaction_file = os.path.join(args.outputdir, 'Reactions.tsv')
    make_tsv(reactions, reaction_file, reaction_header)

    # Write JSON
    reaction_json = open(os.path.join(args.outputdir, 'Reactions.json'), 'w')
    json.dump(reactions, reaction_json, sort_keys=True, indent=4)

    print('Stored reactions in ' + reaction_file + ' and Reactions.json')
