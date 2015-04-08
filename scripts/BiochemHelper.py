
import re

''' Helper methods for working with reactions and compounds in Biochemistry objects. '''

class BiochemHelper:
    
    def __init__(self):
        ''' Initialize object.
        '''

        pass

    def readCompoundsFile(self, path, includeLinenum=True):
        ''' Read the contents of a compounds file.
    
            There is one compound per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.
    
            @param path: Path to compounds file
            @param includeLinenum: When True, include line number in dictionary
            @return List of compound dictionaries.
        '''
    
        # Read the compounds from the specified file.
        compounds = list()
        with open(path, 'r') as handle:
            # The first line has the header with the field names.
            nameList = handle.readline().strip().split('\t')
            fieldNames = dict()
            for index in range(len(nameList)):
                fieldNames[nameList[index]] = index
            required = { 'id', 'name', 'abbreviation', 'formula', 'charge', 'isCofactor' }
            for req in required:
                if req not in fieldNames:
                    print 'WARNING: Required field %s is missing from header' %(req)
                    return None
            
            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.strip().split('\t')
                if len(fields) < len(fieldNames):
                    print 'WARNING: Compound on line %d is missing one or more fields, %s' %(linenum, fields)
                    continue
                cpd = dict()
                cpd['id'] = fields[fieldNames['id']]
                cpd['name'] = fields[fieldNames['name']]
                cpd['formula'] = fields[fieldNames['formula']]
                cpd['charge'] = int(fields[fieldNames['charge']])
                if includeLinenum:
                    cpd['linenum'] = linenum
                compounds.append(cpd)
    
        return compounds
    
    def readReactionsFile(self, path):
        ''' Read the contents of a reactions file.
    
            There is one reaction per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.
    
            @param path: Path to reactions file
            @return List of reaction dictionaries.
        '''
    
        # Read the reactions from the specified file.
        reactions = list()
        with open(path, 'r') as handle:
            # The first line has the header with the field names.
            fieldNames = handle.readline().strip().split('\t')
            # @todo Validate required fields are in the header
            
            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.strip('\n ').split('\t')
                if len(fields) < len(fieldNames):
                    print 'WARNING: Reaction on line %d is missing one or more fields, %s' %(linenum, fields)
                    continue
                rxn = dict()
                rxn['id'] = fields[0]
                rxn['name'] = fields[1]
                rxn['equation'] = fields[2]
                rxn['definition'] = fields[3]
                rxn['status'] = fields[4]
                rxn['linenum'] = linenum
                reactions.append(rxn)
    
        return reactions
    
    def buildDictFromListOfObjects(self, objectList, key='id'):
        ''' Build a dictionary with the specified key from a list of objects.
    
            The value of each element in the returned dictionary is the complete
            object.  This copies the objects into a separate data structure.
    
            @param objectList: List of objects
            @return Dictionary mapping key to object
        '''
    
        objectDict = dict()
        for index in range(len(objectList)):
            obj = objectList[index]
            if obj is not None:
                objectDict[obj[key]] = obj
        return objectDict

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
        if stoichString.find('(') >= 0:
            compound['stoich'] = re.sub(r'\(([0-9]*\.?[0-9]+)\) .+', r'\1', stoichString)
            prefix = stoichString.find(')')
            id = stoichString[prefix+2:] # Maybe this should be converted to a number?
        else:
            compound['stoich'] = 1
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
        if stoichString.find('(') >= 0:
            compound['stoich'] = re.sub(r'\(([0-9]*\.?[0-9]+)\).+', r'\1', stoichString) # MBM remove space?
            prefix = stoichString.find(')')
            id = stoichString[prefix+1:] # Maybe this should be converted to a number?
        else:
            compound['stoich'] = 1
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
        leftdirectional = delimiter+'<='+delimiter
        rightdirectional = delimiter+'=>'+delimiter
        separator = delimiter+'+'+delimiter
    
        # Find the special string that separates reactants and products.
        reactants = list()
        products = list()
        if equation.find(rightdirectional) >= 0:
            direction = '>'
            parts = equation.split(rightdirectional)
            if parts[0]:
                reactants = parts[0].split(separator)
            if parts[1]:
                products = parts[1].split(separator)
        elif equation.find(leftdirectional) >= 0:
            direction = '<'
            parts = equation.split(leftdirectional)
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

            A transport reaction is defined as a reaction in which there is a
            compound on both sides of the reaction and the compound's compartment
            on the reactant side is different from the compound's compartment on
            the product side.

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
            print 'What the hell? '+equation
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
        # the compounds.
        reactantCompounds = dict()
        for cpd in reactants:
            if byName is True:
                compound = self.parseCompoundNameStoich(cpd)
            else:
                compound = self.parseCompoundIdStoich(cpd)
            reactantCompounds[compound['compound']] = compound

        productCompounds = dict()
        for cpd in products:
            if byName is True:
                compound = self.parseCompoundNameStoich(cpd)
            else:
                compound = self.parseCompoundIdStoich(cpd)
            productCompounds[compound['compound']] = compound

        # Run through the list of reactant compounds and see if the compound is
        # also in the list of product compounds.  If so and the compartment is 
        # different then this is a transport reaction.
        transportReaction = False
        transporterList = list()
        for compound in reactantCompounds:
            if compound in productCompounds:
                if reactantCompounds[compound]['compartmentId'] != productCompounds[compound]['compartmentId']:
                    transportReaction = True
                    transporterList.append(reactantCompounds[compound]) # @todo Not sure about what to put in the list

        return transportReaction

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
