#! /usr/bin/env python

import argparse
import os

desc1 = '''
NAME
      Build_Complex_File -- build a complexes file from source files

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Build a complexes file for model templates from source files. The source
      files for ModelSEED and PlantSEED are created by running the
      Print_ModelTemplate_from_Workspaces.pl script.  The source file for KEGG
      is created by running Print_KEGG_Complex_Role.py.
      
      The IDs for complexes and roles are merged into one ID space. ModelSEED
      IDs start at 1, KEGG IDs start at 30000, and PlantSEED IDs start at 50000.
'''

desc3 = '''
EXAMPLES
      Build a complexes file:
      > Build_Complex_File.py Complexes.tsv
      
SEE ALSO
      Build_Role_File.py

AUTHORS
      Mike Mundy 
'''

def readMappingFile(sourceFile, source, rolePrefix, adjust, output):
    ''' Read a Mapping_Complexes file.
    
        @param sourceFile: Path to source Mapping_Complexes file
        @param source: Type of source file (ModelSEED or PlantSEED)
        @param rolePrefix: Prefix string in role references
        @param adjust: Value for adjusting IDs
        @param output: List of output lines for destination file
        @return Nothing
    '''
    
    with open(sourceFile, 'r') as s:
        header = s.readline() # Ignore header line
        for line in s:
            sfields = line.strip().split('\t')
            dfields = list()
            # Need special handling for PlantSEED source which has empty field 1
            if source is 'PlantSEED':
                cid = int(sfields[0].split('.')[1]) + adjust
            else:
                cid = int(sfields[1].strip('cpx')) + adjust
            dfields.append('cpx%05d' %(cid)) # First field is complex ID
            dfields.append(sfields[0]) # Second field is complex name
            dfields.append(source) # Third field is source type
            dfields.append('null') # Fourth field is reference (never specified)
            dfields.append('1.0') # Fifth field is confidence value
            if len(sfields) > 2 and sfields[2] != '': # Sixth field is list of roles
                sourceRoles = sfields[2].split('|')
                destRoles = list()
                for index in range(len(sourceRoles)):
                    sourceRole = sourceRoles[index].split(';')
                    # Assume that roles in source are in a consistent format.
                    values = list()
                    if sourceRole[0].startswith('role_ref:'):
                        # Feature ids
                        fid = int(sourceRole[0].split('/')[-1].strip(rolePrefix)) + adjust
                        values.append('ftr%05d' %(fid))
                    else:
                        print(('Bad role '+sourceRole))
                    if sourceRole[2].startswith('type:'):
                        values.append(sourceRole[2].split(':')[-1])
                    else:
                        print(('Bad role '+sourceRole))
                    if sourceRole[1].startswith('optionalRole:'):
                        values.append(sourceRole[1].split(':')[-1])
                    else:
                        print(('Bad role '+sourceRole))
                    if sourceRole[3].startswith('triggering'):
                        values.append(sourceRole[3].split(':')[-1])
                    else:
                        print(('Bad role '+sourceRole))
                    destRoles.append(';'.join(values))
                dfields.append('|'.join(destRoles))
            else:
                dfields.append('null')
            output.append( {'id': dfields[0], 'line': '\t'.join(dfields) } ) # Makes it easy to sort by ID
    
    return

if __name__ == '__main__':
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Build_Complexes', epilog=desc3)
    parser.add_argument('dest', help='path to merged complexes file', action='store')
    parser.add_argument('--mappingdir', help='path to directory containing mapping source files', action='store', default='../Mappings')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Add each line of the destination complexes file to a list.    
    output = list()
    
    # Get the complexes from the ModelSEED and PlantSEED mapping objects.
    readMappingFile(os.path.join(args.mappingdir, 'default-mapping/Mapping_Complexes.txt'), 'ModelSEED', 'msfr.', 0, output)
    readMappingFile(os.path.join(args.mappingdir, 'PlantSEED_Mapping/Mapping_Complexes.txt'), 'PlantSEED', 'psfr.', 50000, output)
    
    # Get the complexes from the KEGG complex and roles file.         
    source = os.path.join(args.mappingdir, 'ComplexRoles.kegg.tsv')     
    with open(source, 'r') as s:
        header = s.readline()
        for line in s:
            sfields = line.strip().split('\t')
            dfields = list()
            cid = int(sfields[0].split('.')[1]) + int(20000)
            dfields.append('cpx%05d' %(cid)) # First field is complex ID
            dfields.append(sfields[1]) # Second field is complex name
            dfields.append('KEGG') # Third field is source type
            dfields.append('null') # Fourth field is reference (never specified)
            dfields.append('1.0') # Fifth field is confidence value
            values = list() # Sixth field is list of roles
            fid = int(sfields[4].split('.')[1]) - int(20000)
            values.append('ftr%05d' %(fid))
            values.append('triggering')
            values.append(sfields[12])
            values.append(sfields[11])
            dfields.append(';'.join(values))
            output.append( {'id': dfields[0], 'line': '\t'.join(dfields) } ) # Makes it easy to sort by ID

    # Sort the output lines by complex ID.         
    output.sort(key=lambda k: k['id'])
    
    # Write the destination file.
    with open(args.dest, 'w') as d:
        d.write('%s\n' %('\t'.join(['id', 'name', 'source', 'reference', 'confidence', 'roles'])))
        for index in range(len(output)):
            d.write('%s\n' %(output[index]['line']))

    exit(0)