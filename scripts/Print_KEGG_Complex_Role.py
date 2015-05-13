#! /usr/bin/env python

import argparse
import re
from biokbase.probabilistic_annotation.kegg.KEGGEnzymeDatabase import KEGGEnzymeDatabase

desc1 = '''
NAME
      Print_KEGG_Complex_Role -- 

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
'''

desc3 = '''
EXAMPLES
      
SEE ALSO
      Print_Complex_Roles.pl

AUTHORS
      Mike Mundy 
'''

def nameToSearchName(name):
    
    searchName = name.lower()
    searchName = re.sub(r'[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+', r'', searchName) # Remove EC number
    searchName = re.sub(r'\s', r'', searchName) # Remove whitespace
    searchName = re.sub(r'\#.*$', r'', searchName) # Remove special characters
    return searchName

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Print_KEGG_Complex_Role', epilog=desc3)
    parser.add_argument('--enzymedb', help='path to file with KEGG enzyme database', action='store', dest='enzymeDB', default='enzyme.db')
    parser.add_argument('--directory', help='path to directory for storing complex role file', action='store', dest='directory', default='.')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Load the enzyme database.    
    enzymeDB = KEGGEnzymeDatabase(args.enzymeDB)
    enzymeDB.load()
    
    # Assign a complex ID to every enzyme in the enzyme database.
    nextComplexId = 10000
    enzymeComplexDict = dict()
    for key in sorted(enzymeDB.enzymes):
        enzymeComplexDict[key] = 'cpx.%d' %(nextComplexId)
        nextComplexId += 1

    # Build a dictionary keyed by complex ID with a list of roles.
    # We call the enzyme the complex and the roles are the names of the enzyme.
    complexToRoles = dict()
    for key in enzymeComplexDict:
        enzyme = enzymeDB.get(key)
        if len(enzyme.name) > 0 and not enzyme.name[0].startswith('Transferred to') and not enzyme.name[0].startswith('Deleted entry'):
            complexToRoles[enzymeComplexDict[key]] = list()
            for nindex in range(len(enzyme.name)):
                complexToRoles[enzymeComplexDict[key]].append('%s (EC %s)' %(enzyme.name[nindex], key))
#         else:
#             complexToRoles[enzymeComplexDict[key]].append('Unknown (EC %s)' %(key))

    # Build the list of field names in the complex role file.
    fieldNames = [ 
        'complex_id',
        'complex_name',
        'complex_source',
        'complex_type',
        'role_id',
        'role_name',
        'role_type',
        'role_source',
        'role_aliases',
        'role_exemplar',
        'type',
        'triggering',
        'optional'
    ]
    
    # Generate the complex role file.
    nextRoleId = 50000
    with open(args.directory+'/ComplexRoles.kegg.tsv', 'w') as handle:
        handle.write('\t'.join(fieldNames)+'\n')
        for complex in complexToRoles:
            data = list()
            data.append(complex)
            data.append(complex)
            data.append('KEGG')
            data.append('KEGG_role_complex')
            data.append('fr.%d' %(nextRoleId))
            nextRoleId += 1
            data.append(complexToRoles[complex][0]) # Use the first name from KEGG as the role name
            data.append('KEGG_role')
            data.append('KEGG')
            aliases = list()
            aliases.append('searchname:'+nameToSearchName(complexToRoles[complex][0]))
            for index in range(1, len(complexToRoles[complex])):
                aliases.append('kegg:'+complexToRoles[complex][index])
            data.append(';'.join(aliases))
            data.append('null')
            data.append('role_mapping')
            data.append('1')
            data.append('0')
            handle.write('\t'.join(data)+'\n')
