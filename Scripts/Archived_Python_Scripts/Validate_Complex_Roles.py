#! /usr/bin/env python

import argparse
import re
from BiochemHelper import BiochemHelper

desc1 = '''
NAME
      Validate_Complex_Roles -- validate a complex role mapping file

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Validate a complex role mapping file by checking for duplicate IDs and
      names.  The cpxrolefile argument specifies the path to a complex role
      file which describes the mapping from complexes to functional roles in a
      biochemistry.  The first line of the file is a header with the field
      names.  There is one mapping per line with fields separated by tabs.

      When the --show-details optional argument is specified, details on all
      problems are displayed.  When the other --show optional arguments are
      specified, details on the corresponding type of problem are displayed.
'''

desc3 = '''
EXAMPLES
      Show summary data about a complex role mapping file:
      > Validate_Complex_Roles.py ComplexRoles.tsv
      
      Show details on complex role mappings that have problems:
      > Validate_Complex_Roles.py --show-details ComplexRoles.tsv

SEE ALSO
      Validate_Reactions
      Validate_Compounds

AUTHORS
      Mike Mundy 
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Validate_Complex_Roles', epilog=desc3)
    parser.add_argument('cpxrolefile', help='path to complex role mapping file', action='store')
    parser.add_argument('--show-details', help='show details on all problems', action='store_true', dest='showDetails', default=False)
    parser.add_argument('--show-dup-ids', help='show details on duplicate IDs', action='store_true', dest='showDupIds', default=False)
    parser.add_argument('--show-bad-ids', help='show details on bad IDs', action='store_true', dest='showBadIds', default=False)
    parser.add_argument('--show-dup-names', help='show details on duplicate names', action='store_true', dest='showDupNames', default=False)
    parser.add_argument('--show-bad-names', help='show details on bad names', action='store_true', dest='showBadNames', default=False)
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

    # Read the complex role mapping from the specified file.
    print('Complex role mapping file: %s' % args.cpxrolefile)
    helper = BiochemHelper()
    complexRoles = helper.readComplexRolesFile(args.cpxrolefile)
    if complexRoles is None:
        print('Error reading complex role mapping file')
        exit(1)
    print('Number of mappings: %d' % len(complexRoles))

    # Check for duplicates, missing and invalid values.
    complexIdDict = dict()
    duplicateComplexId = 0
    badComplexIdChars = list()
    roleIdDict = dict()
    duplicateRoleId = 0
    badRoleIdChars = list()
    complexNameDict = dict()
    duplicateComplexName = 0
    badComplexNameChars = list()
    roleNameDict = dict()
    duplicateRoleName = 0
    badRoleNameChars = list()
    for index in range(len(complexRoles)):
        cpxrole = complexRoles[index]
        
        # Check for duplicate IDs.
        if cpxrole['complex_id'] in complexIdDict:
            if len(complexIdDict[cpxrole['complex_id']]) == 1:
                duplicateComplexId += 1
            complexIdDict[cpxrole['complex_id']].append(index)
        else:
            complexIdDict[cpxrole['complex_id']] = [index]
        if cpxrole['role_id'] in complexIdDict:
            if len(complexIdDict[cpxrole['role_id']]) == 1:
                duplicateComplexId += 1
            complexIdDict[cpxrole['role_id']].append(index)
        else:
            complexIdDict[cpxrole['role_id']] = [index]

        # Check for invalid characters in the ID.
        match = re.search(r'cpx.\d+', cpxrole['complex_id'])
        if match is None:
            badComplexIdChars.append(index)
        match = re.search(r'fr.\d+', cpxrole['role_id'])
        if match is None:
            badRoleIdChars.append(index)

        # Check for duplicate names.
        if cpxrole['complex_name'] in complexNameDict:
            if len(complexNameDict[cpxrole['complex_name']]) == 1:
                duplicateComplexName += 1
            complexNameDict[cpxrole['complex_name']].append(index)
        else:
            complexNameDict[cpxrole['complex_name']] = [index]
        if cpxrole['role_name'] in roleNameDict:
            if len(roleNameDict[cpxrole['role_name']]) == 1:
                duplicateRoleName += 1
            roleNameDict[cpxrole['role_name']].append(index)
        else:
            roleNameDict[cpxrole['role_name']] = [index]
 
        # Check for invalid characters in the name.
        try:
            cpxrole['complex_name'].encode('ascii')
        except UnicodeEncodeError:
            badComplexNameChars.append(index)
        try:
            cpxrole['role_name'].encode('ascii')
        except UnicodeEncodeError:
            badRoleNameChars.append(index)

    # Print summary data.
    print('Number of mappings with duplicate complex IDs: %d' % duplicateComplexId)
    print('Number of mappings with role complex IDs: %d' % duplicateRoleId)
    print('Number of mappings with bad characters in complex ID: %d' % len(badComplexIdChars))
    print('Number of mappings with bad characters in role ID: %d' % len(badComplexIdChars))
    print('Number of mappings with duplicate complex names: %d' % duplicateComplexName)
    print('Number of mappings with duplicate role names: %d' % duplicateRoleName)
    print('Number of mappings with bad characters in complex name: %d' % len(badComplexNameChars))
    print('Number of mappings with bad characters in role name: %d' % len(badRoleNameChars))
    print()

    # Print details if requested.
    if args.showDupIds:
        for id in complexIdDict:
            if len(complexIdDict[id]) > 1:
                print('Duplicate complex ID: %s' % id)
                for dup in complexIdDict[id]:
                    print('Line %05d: %s' % (complexRoles[dup]['linenum'], complexRoles[dup]))
                print()
        for id in roleIdDict:
            if len(roleIdDict[id]) > 1:
                print('Duplicate role ID: %s' % id)
                for dup in roleIdDict[id]:
                    print('Line %05d: %s' % (complexRoles[dup]['linenum'], complexRoles[dup]))
                print()
    if args.showBadIds:
        if len(badComplexIdChars) > 0:
            print('Complex IDs with bad characters:')
            for index in range(len(badComplexIdChars)):
                print('Line %05d: %s' % (complexRoles[badComplexIdChars[index]]['linenum'], complexRoles[badComplexIdChars[index]]))
            print()
        if len(badRoleIdChars) > 0:
            print('Role IDs with bad characters:')
            for index in range(len(badRoleIdChars)):
                print('Line %05d: %s' % (complexRoles[badRoleIdChars[index]]['linenum'], complexRoles[badRoleIdChars[index]]))
            print()
    if args.showDupNames:
        for name in complexNameDict:
            if len(complexNameDict[name]) > 1:
                print('Duplicate complex name: %s' % name)
                for dup in complexNameDict[name]:
                    print('Line %05d: %s' % (complexRoles[dup]['linenum'], complexRoles[dup]))
                print()
        for name in roleNameDict:
            if len(roleNameDict[name]) > 1:
                print('Duplicate role name: %s' % name)
                for dup in roleNameDict[name]:
                    print('Line %05d: %s' % (complexRoles[dup]['linenum'], complexRoles[dup]))
                print()
    if args.showBadNames:
        if len(badComplexNameChars) > 0:
            print('Complex names with bad characters:')
            for index in range(len(badComplexNameChars)):
                print('Line %05d: %s' % (complexRoles[badComplexNameChars[index]]['linenum'], complexRoles[badComplexNameChars[index]]))
            print()
        if len(badRoleNameChars) > 0:
            print('Role names with bad characters:')
            for index in range(len(badRoleNameChars)):
                print('Line %05d: %s' % (complexRoles[badRoleNameChars[index]]['linenum'], complexRoles[badRoleNameChars[index]]))
            print()

    if any([duplicateComplexId, duplicateRoleId, badComplexIdChars,
            duplicateComplexName, duplicateRoleName, badComplexNameChars,
            badRoleNameChars]):
        exit(0)
