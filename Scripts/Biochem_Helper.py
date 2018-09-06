import re
import os
import sys
from Scripts.Base_Helper import BaseHelper

''' Helper methods for working with Biochemistry source files. '''

class BiochemHelper(BaseHelper):
    
    def __init__(self):
        ''' Initialize object.
        '''

        pass

    def readCompoundsFile(self, path, includeLinenum=True, noFormat=False):
        ''' Read the contents of a compounds file.
    
            There is one compound per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.
            
            The returned compound dictionary uses the field names as keys.  When
            noFormat is True, the field values are exactly as read from the file.
            Otherwise, the fields with a null value are not included in the
            dictionary and numeric values are converted to numbers.
    
            @param path: Path to compounds file
            @param includeLinenum: When True, include line number in dictionary
            @param noFormat: When True, values in compound dictionary are not formatted
            @return List of compound dictionaries.
        '''
    
        # The following fields are required in a compounds file.
        required = { 'id', 'abbreviation', 'name', 'formula', 'mass', 'source',
                      'smiles', 'charge', 'is_core', 'is_obsolete', 'linked_compound',
                      'is_cofactor', 'deltag', 'deltagerr', 'pka', 'pkb',
                      'abstract_compound', 'comprised_of', 'aliases' }

        # Read the compounds from the specified file.
        compounds = list()
        with open(path, 'rU') as handle:
            # The first line has the header with the field names.
            nameList = handle.readline().strip().split('\t')
            fieldNames = self.validateHeader(nameList, required)
            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.rstrip('\n').split('\t',len(fieldNames))
                if len(fields) < len(fieldNames):
                    print('WARNING: Compound on line %d is missing one or more fields, %s' %(linenum, fields))
                    continue
                cpd = dict()
                if noFormat:
                    for index in range(len(nameList)):
                        if(fields[index]=='null'):
                            fields[index] = None
                        if(nameList[index] == 'smiles'):
                            cpd['structure'] = fields[index]
                        elif(nameList[index] == 'inchikey'):
                            pass
                        else:
                            cpd[nameList[index]] = fields[index]
                else:
                    cpd['id'] = fields[fieldNames['id']]
                    cpd['abbreviation'] = fields[fieldNames['abbreviation']]
                    cpd['name'] = fields[fieldNames['name']]
                    cpd['formula'] = fields[fieldNames['formula']]
                    if fields[fieldNames['mass']] != 'null':
                        cpd['mass'] = fields[fieldNames['mass']]
                    else:
                        cpd['mass'] = float(10000000)
                    cpd['source'] = fields[fieldNames['source']]
                    cpd['structure'] = fields[fieldNames['smiles']]
                    if fields[fieldNames['charge']] != 'null':
                        cpd['charge'] = float(fields[fieldNames['charge']])
                    else:
                        cpd['charge'] = float(0)
                    cpd['is_core'] = int(fields[fieldNames['is_core']])
                    cpd['is_obsolete'] = int(fields[fieldNames['is_obsolete']])
                    if fields[fieldNames['linked_compound']] != 'null':
                        cpd['linked_compound'] = fields[fieldNames['linked_compound']]
                    cpd['is_cofactor'] = int(fields[fieldNames['is_cofactor']])
                    if fields[fieldNames['deltag']] != 'null':
                        cpd['deltag'] = float(fields[fieldNames['deltag']])
                    else:
                        cpd['deltag'] = float(10000000)
                    if fields[fieldNames['deltagerr']] != 'null':
                        cpd['deltagerr'] = float(fields[fieldNames['deltagerr']])
                    else:
                        cpd['deltagerr'] = float(10000000)
                    cpd['pka'] = fields[fieldNames['pka']]
                    cpd['pkb'] = fields[fieldNames['pkb']]
                    if fields[fieldNames['abstract_compound']] != 'null':
                        cpd['abstract_compound'] = fields[fieldNames['abstract_compound']]
                    if fields[fieldNames['comprised_of']] != 'null':
                        cpd['comprised_of'] = fields[fieldNames['comprised_of']]
                    cpd['aliases'] = list()
                    if fields[fieldNames['aliases']] != 'null':
                        # This actually needs to parse the aliases format.
                        cpd['aliases'] = fields[fieldNames['aliases']]
                if includeLinenum:
                    cpd['linenum'] = linenum
                compounds.append(cpd)
    
        return compounds
    
    def readReactionsFile(self, path, includeLinenum=True, noFormat=False):
        ''' Read the contents of a reactions file.
    
            There is one reaction per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.
    
            The returned reaction dictionary uses the field names as keys.  When
            noFormat is True, the field values are exactly as read from the file.
            Otherwise, the fields with a null value are not included in the
            dictionary and numeric values are converted to numbers.
    
            @param path: Path to reactions file
            @param includeLinenum: When True, include line number in dictionary
            @param noFormat: When True, values in reaction dictionary are not formatted
            @return List of reaction dictionaries.
        '''

        # The following fields are required in a reactions file.
        required = { 'id', 'abbreviation', 'name', 'code', 'stoichiometry', 'is_transport', 
                     'equation', 'definition', 'reversibility', 'direction', 'abstract_reaction',
                     'pathways', 'aliases', 'ec_numbers', 'deltag', 'deltagerr', 'compound_ids',
                     'status' }

        # Read the reactions from the specified file.
        reactions = list()
        with open(path, 'r') as handle:
            # The first line has the header with the field names.
            nameList = handle.readline().strip().split('\t')
            fieldNames = self.validateHeader(nameList, required)

            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.strip('\n ').split('\t')
                if len(fields) < len(fieldNames):
                    print('WARNING: Reaction on line %d is missing one or more fields, %s' %(linenum, fields))
                    continue
                rxn = dict()
                if noFormat:
                    for index in range(len(nameList)):
                        if(fields[index]=='null'):
                            fields[index] = None
                        if(nameList[index] == 'notes'):
                            pass
                        else:
                            rxn[nameList[index]] = fields[index]
                else:
                    rxn['id'] = fields[fieldNames['id']]
                    rxn['abbreviation'] = fields[fieldNames['abbreviation']]
                    rxn['name'] = fields[fieldNames['name']]
                    rxn['code'] = fields[fieldNames['code']]
                    rxn['stoichiometry'] = fields[fieldNames['stoichiometry']]
                    rxn['is_transport'] = int(fields[fieldNames['is_transport']])
                    rxn['equation'] = fields[fieldNames['equation']]
                    rxn['definition'] = fields[fieldNames['definition']]
                    rxn['reversibility'] = fields[fieldNames['reversibility']]
                    rxn['direction'] = fields[fieldNames['direction']]
                    rxn['abstract_reaction'] = fields[fieldNames['abstract_reaction']]
                    rxn['pathways'] = fields[fieldNames['pathways']]
                    rxn['aliases'] = fields[fieldNames['aliases']]
                    rxn['ec_numbers'] = fields[fieldNames['ec_numbers']]
                    if fields[fieldNames['deltag']] != 'null':
                        rxn['deltag'] = float(fields[fieldNames['deltag']])
                    else:
                        rxn['deltag'] = float(10000000)
                    if fields[fieldNames['deltagerr']] != 'null':
                        rxn['deltagerr'] = float(fields[fieldNames['deltagerr']])
                    else:
                        rxn['deltagerr'] = float(10000000)
                    rxn['compound_ids'] = fields[fieldNames['compound_ids']]
                    rxn['status'] = fields[fieldNames['status']]
                    if 'is_obsolete' in fieldNames:
                        rxn['is_obsolete'] = int(fields[fieldNames['is_obsolete']])
                    if 'linked_reaction' in fieldNames:
                        if fields[fieldNames['linked_reaction']] != 'null':
                            rxn['linked_reaction'] = fields[fieldNames['linked_reaction']]
                if includeLinenum:
                    rxn['linenum'] = linenum
                reactions.append(rxn)
    
        return reactions
    
    def readCompartmentsFile(self, path, includeLinenum=True):
        ''' Read the contents of a compartments file.

            There is one compartment per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.

            @param path: Path to compartments file
            @param includeLinenum: When True, include line number in dictionary
            @return List of compartment dictionaries.
            
            @note This method should become obsolete when compartments are moved to Model Template
        '''

        # The following fields are required in a compartments file.
        required = { 'id', 'name', 'hierarchy' }

        # Read the compartments from the specified file.
        compartments = list()
        with open(path, 'r') as handle:
            # The first line has the header with the field names.
            nameList = handle.readline().strip().split('\t')
            fieldNames = self.validateHeader(nameList, required)

            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.strip('\n ').split('\t')
                if len(fields) < len(fieldNames):
                    print('WARNING: Compartment on line %d is missing one or more fields, %s' %(linenum, fields))
                    continue
                cmp = dict()
                cmp['id'] = fields[fieldNames['id']]
                cmp['name'] = fields[fieldNames['name']]
                cmp['hierarchy'] = int(fields[fieldNames['hierarchy']])
                if includeLinenum:
                    cmp['linenum'] = linenum
                compartments.append(cmp)

        return compartments

    def readComplexRolesFile(self, path, includeLinenum=True):
        ''' Read the contents of a complex role mapping file.

            There is one mapping per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.

            @param path: Path to complex role mapping file
            @param includeLinenum: When True, include line number in dictionary
            @return List of complex role mapping dictionaries.
        '''

        # The following fields are required in a complex roles files.
        required = { 'complex_id', 'complex_name', 'complex_source', 'complex_type', 'role_id',
                     'role_name', 'role_type', 'role_source', 'role_aliases', 'role_exemplar',
                     'type', 'triggering', 'optional' }

        # Read the complex role mappings from the specified file.
        complexRoles = list()
        with open(path, 'r') as handle:
            # The first line has the header with the field names.
            nameList = handle.readline().strip().split('\t')
            fieldNames = self.validateHeader(nameList, required)

            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.strip('\n ').split('\t')
                if len(fields) < len(fieldNames):
                    print('WARNING: Complex role mapping on line %d is missing one or more fields, %s' %(linenum, fields))
                    continue
                cpxrole = dict()
                cpxrole['complex_id'] = fields[fieldNames['complex_id']]
                cpxrole['complex_name'] = fields[fieldNames['complex_name']]
                cpxrole['complex_source'] = fields[fieldNames['complex_source']]
                cpxrole['complex_type'] = fields[fieldNames['complex_type']]
                cpxrole['role_id'] = fields[fieldNames['role_id']]
                cpxrole['role_name'] = fields[fieldNames['role_name']]
                cpxrole['role_type'] = fields[fieldNames['role_type']]
                cpxrole['role_source'] = fields[fieldNames['role_source']]
                cpxrole['role_aliases'] = fields[fieldNames['role_aliases']]
                cpxrole['role_exemplar'] = fields[fieldNames['role_exemplar']]
                cpxrole['type'] = fields[fieldNames['type']]
                cpxrole['triggering'] = int(fields[fieldNames['triggering']])
                cpxrole['optional'] = int(fields[fieldNames['optional']])
                if includeLinenum:
                    cpxrole['linenum'] = linenum
                complexRoles.append(cpxrole)

        return complexRoles

    def readAliasFiles(self, aliasDir):
        ''' Read the contents of all of the alias files.
        
            There is one alias mapping per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.
            
            The returned dictionaries use the ModelSEED ID as the key and each source has
            a list of values that are aliases for the compound or reaction.
        
            @param aliasDir: Path to directory containing alias files.
            @return Dictionary of compound aliases, dictionary of reaction aliases.
        '''

        compoundAliases = dict() # First level is compound ID, second level is source with a list of values.
        reactionAliases = dict()
    
        # Read all of the alias files.  Each alias file has a header on the first line. Assumption
        # is that first column is the ID in the alternate system, second column is the ID
        # in ModelSEED, and third column is the ID in PlantSEED.  Multiple values in the second
        # and third columns are separated by "|".  Second or third column can be empty if there
        # is no match in ModelSEED or PlantSEED.  Aliases can be either for compounds or
        # reactions as indicated by the prefix on the ID.
        aliasFiles = os.listdir(aliasDir)
        for index in range(len(aliasFiles)):
            (source, ext) = os.path.splitext(aliasFiles[index]) 
            if ext != '.aliases':
                continue
            print('Processing aliases in '+aliasFiles[index])
            source = source.replace('_', ' ') # Replace the underscores used in file names
            with open(os.path.join(aliasDir, aliasFiles[index]), 'r') as handle:
                header = handle.readline().strip().split('\t')
                for line in handle:
                    fields = line.strip().split('\t')
                    # Need to split on | Need to separate in four checks
                    if fields[1].startswith('cpd'):
                        idList = fields[1].split('|')
                        for index in range(len(idList)):
                            if idList[index] not in compoundAliases:
                                compoundAliases[idList[index]] = dict()
                            if source not in compoundAliases[idList[index]]:
                                compoundAliases[idList[index]][source] = list()
                            compoundAliases[idList[index]][source].append(fields[0])
                    if len(fields) > 2 and fields[2].startswith('cpd'):
                        idList = fields[2].split('|')
                        for index in range(len(idList)):
                            if idList[index] not in compoundAliases:
                                compoundAliases[idList[index]] = dict()
                            if source not in compoundAliases[idList[index]]:
                                compoundAliases[idList[index]][source] = list()
                            compoundAliases[idList[index]][source].append(fields[0])
                    if fields[1].startswith('rxn'):
                        idList = fields[1].split('|')
                        for index in range(len(idList)):
                            if idList[index] not in reactionAliases:
                                reactionAliases[idList[index]] = dict()
                            if source not in reactionAliases[idList[index]]:
                                reactionAliases[idList[index]][source] = list()
                            reactionAliases[idList[index]][source].append(fields[0])
                    if len(fields) > 2 and fields[2].startswith('rxn'):
                        idList = fields[2].split('|')
                        for index in range(len(idList)):
                            if idList[index] not in reactionAliases:
                                reactionAliases[idList[index]] = dict()
                            if source not in reactionAliases[idList[index]]:
                                reactionAliases[idList[index]][source] = list()
                            reactionAliases[idList[index]][source].append(fields[0])

        return compoundAliases, reactionAliases

    def parseCompoundIdStoich(self, stoichString):
        ''' Parse a compound stoichiometry string with IDs into a dictionary.
    
            A stoichiometry string is in the format: (n) cpd[Cm]
            where n is a stoichiometry coefficient, cpd is a compound ID,
            C is a compartment ID and m is compartment index number.  For
            example: (2) cpd00067[c0]
            
            The output is a dictionary with these keys: (1) 'stoich' has the value
            of the coefficient n, (2) 'id' has the value of the compound ID in the
            format cpd_Cm (for example cpd00067_c0), (3) 'compartment' has the
            value of the compartment (for example c0), (4) 'compartmentId' has the
            value of the compartment ID (for example c), (5) 'compartmentIndex' has
            the value of the compartment index (for example 0), and (6) 'compound'
            has the value of the compound ID without the compartment suffix (for
            example cpd00067).
    
            @param stoichString: String representation of compound stoichiometry
            @return Dictionary with components of string as keys
        '''
    
        # Start with an empty dictionary.
        compound = dict()
        
        # Extract the stoichiometry coefficient from inside the parenthesis.
        lparen = stoichString.find('(')
        if lparen >= 0:
            rparen = stoichString.find(')')
            compound['stoich'] = float(stoichString[lparen+1:rparen])
            id = stoichString[rparen+2:]
        else:
            compound['stoich'] = 1.0
            id = stoichString
            
        # Convert the ID from a format with the compartment and index number in
        # square brackets to a format with an underscore.
        parts = id.split('[')
        compound['compound'] = parts[0]
        if len(parts) > 1:
            spos = parts[1].find(']')
            if spos >= 0:
                compound['compartment'] = parts[1][:spos]
            else:
                compound['compartment'] = parts[1]
            compound['compartmentId'] = compound['compartment'][0]
            if len(compound['compartment']) > 1:
                compound['compartmentIndex'] = compound['compartment'][1:]
            else:
                compound['compartmentIndex'] = '0'
        else:
            compound['compartment'] = 'c0' # If not specified put the compound in cytosol
            compound['compartmentId'] = 'c'
            compound['compartmentIndex'] = '0'
        compound['id'] = '%s_%s' %(compound['compound'], compound['compartment'])
        return compound
    
    def parseCompoundNameStoich(self, stoichString):
        ''' Parse a compound stoichiometry string with names into a dictionary.
    
            A stoichiometry string is in the format: (n)cpdname[Cm]
            where n is a stoichiometry coefficient, cpdname is a compound name,
            C is a compartment ID and m is compartment index number.  For
            example: (0.00793965859468043)CoA[c0]
    
            The output is a dictionary with these keys: (1) 'stoich' has the value
            of the coefficient n, (2) 'compound' has the value of the compound name
            (for example CoA), (3) 'compartment' has the value of the compartment (for
            example c0), (4) 'compartmentId' has the value of the compartment ID
            (for example c), (5) 'compartmentIndex' has the value of the compartment
            index (for example 0).
    
            @param stoichString: String representation of compound stoichiometry
            @return Dictionary with components of string as keys
        '''
    
        # Start with an empty dictionary.
        compound = dict()
        
        # Extract the stoichiometry coefficient from inside the parenthesis.
        lparen = stoichString.find('(')
        if lparen >= 0:
            rparen = stoichString.find(')')
            compound['stoich'] = float(stoichString[lparen+1:rparen])
            id = stoichString[rparen+2:]
        else:
            compound['stoich'] = 1.0
            id = stoichString
    
        # Convert the ID from a format with the compartment and index number in
        # square brackets to a format with an underscore.
        parts = id.split('[')
        compound['compound'] = parts[0]
        if len(parts) > 1:
            spos = parts[1].find(']')
            if spos >= 0:
                compound['compartment'] = parts[1][:spos]
            else:
                compound['compartment'] = parts[1]
            compound['compartmentId'] = compound['compartment'][0]
            if len(compound['compartment']) > 1:
                compound['compartmentIndex'] = compound['compartment'][1:]
            else:
                compound['compartmentIndex'] = '0'
        else:
            compound['compartment'] = 'c0' # If not specified put the compound in cytosol
            compound['compartmentId'] = 'c'
            compound['compartmentIndex'] = '0'
        return compound
    
    def parseEquation(self, equation, delimiter=' '):
        ''' Parse a reaction equation string into two lists of compounds.
        
            The compounds in the reaction equation string can be given either as
            compound names or as compound IDs.  The following is an example
            equation using names:
            (1) H+[e0] + (1) L-Tryptophan[e0] => (1) H+[c0] + (1) L-Tryptophan[c0]

            The following is an example using IDs:
            (1) cpd00011[c0] + (1) cpd02255[c0] => (1) cpd00067[c0] + (1) cpd00938[c0]
            
            The string is parsed into a list of reactants and a list of products.
            Each compound in a list is a string that can be parsed into its components
            by parseCompoundIdStoich() or parseCompoundNameStoich().
            
            @params equation: String representing the reaction equation
            @params delimiter: String delimiting components in equation string
            @return List of reactant compound strings, list of product compound strings
        '''
    
        # Build search strings using specified delimiter.
        bidirectional = delimiter+'<=>'+delimiter
        reverse = delimiter+'<='+delimiter
        forward = delimiter+'=>'+delimiter
        separator = delimiter+'+'+delimiter
    
        # Find the special string that separates reactants and products.
        reactants = list()
        products = list()
        if equation.find(forward) >= 0:
            direction = '>'
            parts = equation.split(forward)
            if parts[0]:
                reactants = parts[0].split(separator)
            if parts[1]:
                products = parts[1].split(separator)
        elif equation.find(reverse) >= 0:
            direction = '<'
            parts = equation.split(reverse)
            if parts[1]:
                reactants = parts[1].split(separator)
            if parts[0]:
                products = parts[0].split(separator)
        elif equation.find(bidirectional) >= 0:
            direction = '='
            parts = equation.split(bidirectional)
            if parts[0]:
                reactants = parts[0].split(separator)
            if parts[1]:
                products = parts[1].split(separator)
        else:
            return None, None # Should consider throwing an exception here
    
        return reactants, products

    def isTransportReaction(self, equation):
        ''' Determine if a reaction is a transport reaction.

            A transport reaction is defined as a reaction in which there are
            compounds in more than one compartment. Typically a transport
            reaction is moving the same compound from one compartment to
            another compartment.

            The compounds in the reaction equation string can be given either as
            compound names or as compound IDs.

            @param equation: String representing the reaction equation
            @return True when the reaction is a transport reaction, otherwise False
        '''

        # Parse the reaction's equation to get a list of reactant compounds and
        # a list of product compounds.  The equation is expressed in terms of
        # compound IDs.
        reactants, products = self.parseEquation(equation)
        if len(reactants) == 0 and len(products) == 0:
            print('This reaction has no reactants and no products: '+equation)
            return False

        # Figure out if the compounds are identified by ID or by name.
        if len(reactants) > 0:
            match = re.search(r'cpd\d+', reactants[0])
        else:
            match = re.search(r'cpd\d+', products[0])
        if match is None:
            byName = True
        else:
            byName = False

        # Parse the reactant and product compound lists to get the details on
        # the compounds. Keep track of the compartment IDs for reactants and
        # products.
        compartments = list()
        for cpd in reactants:
            if byName is True:
                compound = self.parseCompoundNameStoich(cpd)
            else:
                compound = self.parseCompoundIdStoich(cpd)
            compartments.append(compound['compartmentId'])

        for cpd in products:
            if byName is True:
                compound = self.parseCompoundNameStoich(cpd)
            else:
                compound = self.parseCompoundIdStoich(cpd)
            compartments.append(compound['comparmentId'])

        # Run through the list of compartment IDs. If all of the compartment IDs
        # are the same then this is not a transport reaction. Otherwise there are
        # multiple compartments and this is a transport reaction.
        for index in range(1,len(compartments)):
            if compartments[0] != compartments[index]:
                return True
        return False

    def isCompoundIdInList(self, compoundId, compoundList):
        ''' Determine if a compound ID is in a list of compounds.

            Use this method to look for a compound in a list returned by parseEquation().

            @param compoundId: ID of compound to search for
            @param compoundList: List of compound strings
        '''

        # Run through the list of compounds and look for a match.
        for cpd in compoundList:
            compound = self.parseCompoundIdStoich(cpd)
            if compound['compound'] == compoundId:
                return True

        return False

    def isCompoundReactant(self, compoundId, equation, direction):
        ''' Determine if a compound is a reactant in a reaction.

            Only reliable way to check is to use compound IDs (otherwise all
            alias names need to be checked).

            @param compoundId: ID of compound to check
            @param equation: String representing the reaction equation
            @param direction: Direction of the reaction
            @return True when the compound is a reactant, otherwise False
        '''

        # Parse the equation to get separate lists of reactants and products.
        reactantList, productList = self.parseEquation(equation)

        # For a bi-directional reaction, also check the product list for the compound.
        # @todo: Should this also check the reversibility field (even though that is at pH 7)?
        if direction == '=':
            if productList is not None and self.isCompoundIdInList(compoundId, productList):
                return True

        # Check the reactant list for the compound.
        if reactantList is None:
            return False
        return self.isCompoundIdInList(compoundId, reactantList)

    def isCompoundProduct(self, compoundId, equation, direction):
        ''' Determine if a compound is a product in a reaction.

            Only reliable way to check is to use compound IDs (otherwise all
            alias names need to be checked).

            @param compoundId: ID of compound to check
            @param equation: String representing the reaction equation
            @param direction: Direction of the reaction
            @return True when the compound is a product, otherwise False
        '''

        # Parse the equation to get separate lists of reactants and products.
        reactantList, productList = self.parseEquation(equation)

        # For a bi-directional reaction, also check the reactant list for the compound.
        # @todo: Should this also check the reversibility field (even though that is at pH 7)?
        if direction == '=':
            if reactantList is not None and self.isCompoundIdInList(compoundId, reactantList):
                return True

        # Check the product list for the compound.
        if productList is None:
            return False
        return self.isCompoundIdInList(compoundId, productList)
