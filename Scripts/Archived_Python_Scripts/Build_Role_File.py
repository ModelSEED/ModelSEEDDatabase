#! /usr/bin/env python

import argparse
import os
import re

desc1 = '''
NAME
      Build_Role_File -- build a roles file from source files

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Build a roles file for model templates from source files. The source
      files for ModelSEED and PlantSEED are created by running the
      Print_ModelTemplate_from_Workspaces.pl script.  The source file for KEGG
      is created by running Print_KEGG_Complex_Role.py.
      
      The IDs for complexes and roles are merged into one ID space. ModelSEED
      IDs start at 1, KEGG IDs start at 30000, and PlantSEED IDs start at 50000.
'''

desc3 = '''
EXAMPLES
      Build a roles file:
      > Build_Role_File.py Roles.tsv
      
SEE ALSO
      Build_Complex_File.py

AUTHORS
      Mike Mundy 
'''

def convertRoleToSearchRole(role):
    ''' Convert a role to a condensed format for better searches.
    
        @param role: Role string
        @return Condensed role string for searches
    '''
    searchRole = role.lower()
    searchRole = re.sub(r'[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+', '', searchRole) # Remove EC number digits
    searchRole = re.sub(r'\s', '', searchRole) # Remove whitespace
    searchRole = re.sub(r'\#.*$', '', searchRole) # Remove comments from the end
    searchRole = re.sub(r'\(ec\)', '', searchRole) # Remove EC number prefix
    return searchRole

def readMappingFile(sourceFile, source, adjust, delimiter, output):
    with open(sourceFile, 'r') as s:
        header = s.readline()
        for line in s:
            sfields = line.strip().split('\t')
            dfields = list()
            rid = int(sfields[0].split('.')[1]) + adjust
            dfields.append('ftr%05d' %(rid)) # First field is role ID
            dfields.append(sfields[1]) # Second field is role name
            dfields.append(source) # Third field is source type
            if len(sfields) > 2 and sfields[2] != 'NONE': # Fourth field is list of features
                dfields.append(sfields[2])
            else:
                dfields.append('null')
            dfields.append(source+':'+sfields[0]+delimiter+'searchname:'+convertRoleToSearchRole(sfields[1])) # Fifth field is list of aliases
            output.append( { 'id': dfields[0], 'line': '\t'.join(dfields) } )  # Makes it easy to sort by ID

    return

if __name__ == '__main__':
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Build_Complexes', epilog=desc3)
    parser.add_argument('dest', help='path to merged complexes file', action='store')
    parser.add_argument('--mappingdir', help='path to directory containing mapping source files', action='store', default='../Mappings')
    parser.add_argument('--delimiter', help='delimiter for separating aliases', action='store', default=';')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()
    
    # Add each line of the destination roles file to a list.    
    output = list()
    
    # Get the roles from the ModelSEED and PlantSEED mapping objects.
    readMappingFile(os.path.join(args.mappingdir, 'default-mapping/Mapping_Roles.txt'), 'ModelSEED', 0, args.delimiter, output)
    readMappingFile(os.path.join(args.mappingdir, 'PlantSEED_Mapping/Mapping_Roles.txt'), 'PlantSEED', 50000, args.delimiter, output)

    # Get the roles from the KEGG complex and roles file.         
    source = os.path.join(args.mappingdir, 'ComplexRoles.kegg.tsv')
    with open(source, 'r') as s:
        header = s.readline()
        for line in s:
            sfields = line.strip().split('\t')
            dfields = list()
            fid = int(sfields[4].split('.')[1]) - 20000
            dfields.append('ftr%05d' %(fid)) # First field is role ID
            dfields.append(sfields[5]) # Second field is role name
            dfields.append('KEGG') # Third field is source type
            dfields.append('null') # Fourth field is list of features
            # @todo Source file has multiple names for roles in aliases field that
            #       are not propogated to the destination file.
            dfields.append('KEGG:'+sfields[4]+args.delimiter+'searchname:'+convertRoleToSearchRole(sfields[5])) # Fifth field is list of aliases
            output.append( { 'id': dfields[0], 'line': '\t'.join(dfields) } ) # Makes it easy to sort by ID

    # Sort the output lines by role ID.         
    output.sort(key=lambda k: k['id'])
    
    # Write the destination file.
    with open(args.dest, 'w') as d:
        d.write('%s\n' %('\t'.join(['id', 'name', 'source', 'features', 'aliases'])))
        for index in range(len(output)):
            d.write('%s\n' %(output[index]['line']))

    exit(0)