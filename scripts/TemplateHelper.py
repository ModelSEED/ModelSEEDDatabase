
import re
from BaseHelper import BaseHelper
from BiochemHelper import BiochemHelper

class CompoundNotFoundError(Exception):
    pass

class ObsoleteCompoundError(Exception):
    pass

class CompoundKeyError(Exception):
    pass

class DuplicateCompartmentError(Exception):
    pass

class CompartmentNotFoundError(Exception):
    pass

class DuplicateRoleError(Exception):
    pass

class RoleNotFoundError(Exception):
    pass

class ReactionNotFoundError(Exception):
    pass

class ReactionFormatError(Exception):
    pass

class DuplicateReactionError(Exception):
    pass

class ObsoleteReactionError(Exception):
    pass

class ComplexNotFoundError(Exception):
    pass

class NoComplexesError(Exception):
    pass

class LinkedCompoundFormatError(Exception):
    pass

''' Helper methods for working with Model Template objects. '''

class TemplateHelper(BaseHelper):
    
    def __init__(self, compoundsPath, reactionsPath):
        ''' Initialize object.
        
            @param compoundsPath: Path to master compounds file
            @param reactionsPath: Path to master reactions file
            @return Nothing
        '''

        # Load the master compounds and reactions from source files.
        self.biochem = BiochemHelper()
        self.masterCompoundsList = self.biochem.readCompoundsFile(compoundsPath, includeLinenum=False)
        self.masterCompounds = self.buildIndexDictFromListOfObjects(self.masterCompoundsList)
        self.masterReactionsList = self.biochem.readReactionsFile(reactionsPath, includeLinenum=False)
        self.masterReactions = self.buildIndexDictFromListOfObjects(self.masterReactionsList)
        
        # Create empty dictionaries for keeping track of items to add to Model Template.
        self.compartments = dict()
        self.biomasses = dict()
        self.compounds = dict()
        self.compCompounds = dict()
        self.roles = dict()
        self.complexes = dict()
        self.reactions = dict()

        return

    def readBiomassesFile(self, biomassPath, compoundsPath, includeLinenum=True, noFormat=False):
        ''' Read the contents of a biomasses and biomass compounds files.
    
            There is one compound per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.
            
            The TemplateBiomass structure uses the field names as keys.  When
            noFormat is True, the field values are exactly as read from the file.
            Otherwise, the fields with a null value are converted to default values
            and numeric values are converted to numbers.
    
            @param biomassPath: Path to biomasses file
            @param compoundsPath: Path to biomass compounds file
            @param includeLinenum: When True, include line number in dictionary
            @param noFormat: When True, values in compound dictionary are not formatted
            @return Nothing
        '''

        # Read the biomass compounds from the specified file.
        compounds = dict()
        with open(compoundsPath, 'r') as handle:
            # The first line has the header with the field names.
            nameList = handle.readline().strip().split('\t')
            required = { 'biomass_id', 'id', 'coefficient', 'coefficient_type', 'class', 
                         'linked_compounds', 'compartment' }
            fieldNames = self.validateHeader(nameList, required)
            
            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.strip().split('\t')
                if len(fields) < len(required):
                    print 'WARNING: Biomass compound on line %d is missing one or more fields, %s' %(linenum, fields)
                    continue
                
                # Create a new TemplateBiomassComponent.
                component = dict()
                if noFormat:
                    for index in range(len(nameList)):
                        component[nameList[index]] = fields[index]
                else:
                    compCompound = self.addCompCompound(fields[fieldNames['id']], fields[fieldNames['compartment']])
                    component['class'] = fields[fieldNames['class']]
                    component['templatecompcompound_ref'] = '~/compcompounds/id/'+compCompound['id']
                    component['coefficient_type'] = fields[fieldNames['coefficient_type']]
                    component['coefficient'] = float(fields[fieldNames['coefficient']])
                    component['linked_compound_refs'] = list()
                    component['link_coefficients'] = list()
                    if fields[fieldNames['linked_compounds']] != 'null':
                        linkedCpds = fields[fieldNames['linked_compounds']].split('|')
                        if len(linkedCpds) < 1:
                            raise LinkedCompoundFormatError('Biomass compound %s on line %d has an invalid linked compound field' %(fields[fieldNames['id']], linenum))
                        for lcIndex in range(len(linkedCpds)):
                            parts = linkedCpds[lcIndex].split(':')
                            if len(parts) != 2:
                                raise LinkedCompoundFormatError('Biomass compound %s on line %d has an invalid linked compound field' %(fields[fieldNames['id']], linenum))                                
                            linkCompCompound = self.addCompCompound(parts[0], fields[fieldNames['compartment']])
                            component['linked_compound_refs'].append('~/compcompounds/id/'+linkCompCompound['id'])
                            component['link_coefficients'].append(float(parts[1]))
                if includeLinenum:
                    component['linenum'] = linenum
                if fields[fieldNames['biomass_id']] not in compounds:
                    compounds[fields[fieldNames['biomass_id']]] = list()    
                compounds[fields[fieldNames['biomass_id']]].append(component)
            
        # Read the biomasses from the specified file.
        with open(biomassPath, 'r') as handle:
            # The first line has the header with the field names.
            nameList = handle.readline().strip().split('\t')
            required = { 'id', 'name', 'type', 'other', 'dna', 'rna', 'protein', 'lipid', 
                         'cellwall', 'cofactor', 'energy' }
            fieldNames = self.validateHeader(nameList, required)
            
            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.strip().split('\t')
                if len(fields) < len(required):
                    print 'WARNING: Biomass on line %d is missing one or more fields, %s' %(linenum, fields)
                    continue
                
                # Create a new TemplateBiomass.
                biomass = dict()
                if noFormat:
                    for index in range(len(nameList)):
                        biomass[nameList[index]] = fields[index]
                else:
                    biomass['id'] = fields[fieldNames['id']]
                    biomass['name'] = fields[fieldNames['name']]
                    biomass['type'] = fields[fieldNames['type']]
                    biomass['other'] = float(fields[fieldNames['other']])
                    biomass['dna'] = float(fields[fieldNames['dna']])
                    biomass['rna'] = float(fields[fieldNames['rna']])
                    biomass['protein'] = float(fields[fieldNames['protein']])
                    biomass['lipid'] = float(fields[fieldNames['lipid']])
                    biomass['cellwall'] = float(fields[fieldNames['cellwall']])
                    biomass['cofactor'] = float(fields[fieldNames['cofactor']])
                    biomass['energy'] = float(fields[fieldNames['energy']])
                    biomass['templateBiomassComponents'] = compounds[biomass['id']]
                if includeLinenum:
                    biomass['linenum'] = linenum
                if biomass['id'] not in self.biomasses:
                    self.biomasses[biomass['id']] = biomass
                else:
                    raise DuplicateBiomassError('Biomass %s on line %d is a duplicate' %(biomass['id'], linenum))

        return
    
    def readCompartmentsFile(self, path, includeLinenum=True, noFormat=False):
        ''' Read the contents of a compartments file.
    
            There is one compartment per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.
            
            The TemplateCompartment structure uses the field names as keys.  When
            noFormat is True, the field values are exactly as read from the file.
            Otherwise, the fields with a null value are converted to default values
            and numeric values are converted to numbers.
    
            @param path: Path to compartments file
            @param includeLinenum: When True, include line number in dictionary
            @param noFormat: When True, values in compound dictionary are not formatted
            @return Nothing
        '''

        # The following fields are required in a compartments file.
        required = { 'index', 'id', 'name', 'hierarchy', 'pH', 'aliases' }

        # Read the compartments from the specified file.
        with open(path, 'r') as handle:
            # The first line has the header with the field names.
            nameList = handle.readline().strip().split('\t')
            fieldNames = self.validateHeader(nameList, required)
            
            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.strip().split('\t')
                if len(fields) < len(required):
                    print 'WARNING: Compartment on line %d is missing one or more fields, %s' %(linenum, fields)
                    continue
                
                # Create a new TemplateCompartment.
                compartment = dict()
                if noFormat:
                    for index in range(len(nameList)):
                        compartment[nameList[index]] = fields[index]
                else:
                    compartment['index'] = fields[fieldNames['index']]
                    compartment['id'] = fields[fieldNames['id']]
                    compartment['name'] = fields[fieldNames['name']]
                    compartment['hierarchy'] = int(fields[fieldNames['hierarchy']])
                    compartment['pH'] = float(fields[fieldNames['pH']])
                    compartment['aliases'] = list()
                    if fields[fieldNames['aliases']] != 'null':
                        self.addToList(fields[fieldNames['aliases']], ';', compartment['aliases'])
                if includeLinenum:
                    compartment['linenum'] = linenum
                if compartment['id'] not in self.compartments:
                    self.compartments[compartment['id']] = compartment
                else:
                    raise DuplicateCompartmentError('Compartment %s on line %d is a duplicate' %(compartment['id'], linenum))

        return

    def readRolesFile(self, path, includeLinenum=True, noFormat=False):
        ''' Read the contents of a roles file.

            There is one role per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.
            
            The TemplateRole structure uses the field names as keys.  When noFormat
            is True, the field values are exactly as read from the file. Otherwise,
            the fields with a null value are converted to default values and numeric
            values are converted to numbers.
    
            @param path: Path to roles file
            @param includeLinenum: When True, include line number in dictionary
            @param noFormat: When True, values in TemplateRole dictionary are not formatted
            @return Nothing
        '''

        # The following fields are required in a reactions file.
        required = { 'id', 'name', 'source', 'features', 'aliases' }

        # Read the reactions from the specified file.
        with open(path, 'r') as handle:
            # The first line has the header with the field names.
            nameList = handle.readline().strip().split('\t')
            fieldNames = self.validateHeader(nameList, required)
            
            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.strip().split('\t')
                if len(fields) < len(required):
                    print 'WARNING: Role on line %d is missing one or more fields, %s' %(linenum, fields)
                    continue
                

                # Create a new TemplateCompartment.
                role = dict()
                if noFormat:
                    for index in range(len(nameList)):
                        role[nameList[index]] = fields[index]
                else:
                    role['id'] = fields[fieldNames['id']]
                    role['name'] = fields[fieldNames['name']]
                    role['source'] = fields[fieldNames['source']]
                    role['features'] = list()
                    if fields[fieldNames['features']] != 'null':
                        self.addToList(fields[fieldNames['features']], ';', role['features'])
                    role['aliases'] = list()
                    if fields[fieldNames['aliases']] != 'null':
                        self.addToList(fields[fieldNames['aliases']], ';', role['aliases'])
#                        role['aliases'] = self.makeAliases(fields[fieldNames['aliases']], '///', ':')
                if includeLinenum:
                    role['linenum'] = linenum
                if role['id'] not in self.roles:
                    self.roles[role['id']] = role
                else:
                    raise DuplicateRoleError('Role %s on line %d is a duplicate' %(role['id'], linenum))
        return

    def readComplexesFile(self, path, includeLinenum=True, noFormat=False):
        ''' Read the contents of a complex roles file.

            There is one complex per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.
            
            The TemplateComplex structure uses the field names as keys.  When noFormat
            is True, the field values are exactly as read from the file. Otherwise, the
            fields with a null value are converted to default values and numeric values
            are converted to numbers.
    
            @param path: Path to complex roles file
            @param includeLinenum: When True, include line number in dictionary
            @return Nothing
        '''

        # The following fields are required in a complexes file.
        required = { 'id', 'name', 'source', 'reference', 'confidence', 'roles' }

        # Read the reactions from the specified file.
        with open(path, 'r') as handle:
            # The first line has the header with the field names.
            nameList = handle.readline().strip().split('\t')
            fieldNames = self.validateHeader(nameList, required)
            
            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.strip().split('\t')
                if len(fields) < len(required):
                    print 'WARNING: Role on line %d is missing one or more fields, %s' %(linenum, fields)
                    continue
                
                # Create a new TemplateComplex if needed. The same complex can be paired with multiple roles.
                complex = dict()
                if noFormat:
                    for index in range(len(nameList)):
                        role[nameList[index]] = fields[index]
                else:
                    complex['id'] = fields[fieldNames['id']]
                    complex['name'] = fields[fieldNames['name']]
                    complex['source'] = fields[fieldNames['source']]
                    complex['reference'] = fields[fieldNames['reference']]
                    complex['confidence'] = float(fields[fieldNames['confidence']])
                    complex['complexroles'] = list()
                    if fields[fieldNames['roles']] != 'null':
                        # Link the role and complex through a TemplateComplexRole.
                        links = fields[fieldNames['roles']].split('|')
                        for index in range(len(links)):
                            values = links[index].split(';')
                            if values[0] not in self.roles:
                                raise RoleNotFoundError('Role %s on line %d not found' %(values[0], linenum))
                            complexRole = dict()
                            complexRole['templaterole_ref'] = '~/roles/id/'+values[0] # Need to validate
                            complexRole['optional'] = int(values[2])
                            complexRole['triggering'] = int(values[3])
                            complex['complexroles'].append(complexRole)
                if includeLinenum:
                    complex['linenum'] = linenum
                if complex['id'] not in self.complexes:
                    self.complexes[complex['id']] = complex
                else:
                    raise DuplicateComplexError('Complex %s on line %d is a duplicate' %(role['id'], linenum))
                
        return
    
    def readReactionsFile(self, path, includeLinenum=True, noFormat=False):
        ''' Read the contents of a reactions file.

            There is one reaction per line in the file with fields separated by tabs.
            The first line of the file is a header with the field names.
            
            The TemplateReaction structure uses the field names as keys.  When
            noFormat is True, the field values are exactly as read from the file.
            Otherwise, the fields with a null value are converted to default values
            and numeric values are converted to numbers.
    
            @param path: Path to reactions file
            @param includeLinenum: When True, include line number in dictionary
            @param noFormat: When True, values in reaction dictionary are not formatted
            @return Nothing
        '''

        # The following fields are required in a reactions file.
        required = { 'id', 'compartment', 'direction', 'gfdir', 'type', 'base_cost',
                     'forward_cost', 'reverse_cost', 'complexes' }

        # Read the reactions from the specified file.
        with open(path, 'r') as handle:
            # The first line has the header with the field names.
            nameList = handle.readline().strip().split('\t')
            fieldNames = self.validateHeader(nameList, required)
            
            linenum = 1
            for line in handle:
                linenum += 1
                fields = line.strip().split('\t')
                if len(fields) < len(required):
                    print 'WARNING: Reaction on line %d is missing one or more fields, %s' %(linenum, fields)
                    continue
                
                # Create a new TemplateReaction.
                reaction = dict()
                if noFormat:
                    for index in range(len(nameList)):
                        reaction[nameList[index]] = fields[index]
                else:
                    # Get the reaction from the master list.
                    reactionId = fields[fieldNames['id']]
                    try:
                        masterReaction = self.masterReactionsList[self.masterReactions[reactionId]]
                    except:
                        raise ReactionNotFoundError('Reaction %s not found in master biochemistry' %(reactionId))
                    
                    # Check the reaction status.
                    #if 'OK' not in masterReaction['status']:
                    #    print 'WARNING: Reaction %s has status %s and was skipped' %(masterReaction['id'], masterReaction['status'])
                    #    continue

                    # Check for obsolete reaction.
                    if masterReaction['is_obsolete']:
                        if masterReaction['linked_reaction'] != 'null':
                            # One of the reactions in the list is not obsolete.
                            linkIds = masterReaction['linked_reaction'].split(';')
                            for index in range(len(linkIds)):
                                try:
                                    linkReaction = self.masterReactionsList[self.masterReactions[linkIds[index]]]
                                    if not linkReaction['is_obsolete']:
                                        print 'NOTICE: Obsolete reaction %s replaced by %s' %(masterReaction['id'], linkReaction['id'])
                                        masterReaction = linkReaction
                                        break
                                except KeyError as e:
                                    raise ObsoleteReactionError('Reaction %s is obsolete and replacement %s not found' %(masterReaction['id'], linkIds[index]))
                            if masterReaction['is_obsolete']:
                                raise ObsoleteReactionError('Reaction %s is obsolete and all replacements are obsolete' %(masterReaction['id']))
                        else:
                            raise ObsoleteReactionError('Reaction %s is obsolete and no replacement is specified' %(masterReaction['id']))
                    
                    # Make sure all of the compartments are valid.
                    compartmentIds = fields[fieldNames['compartment']].split('|')
                    idcomp = compartmentIds[0]
                    for cindex in range(len(compartmentIds)):
                        try:
                            self.compartments[compartmentIds[cindex]]
                        except KeyError as e:
                            raise CompartmentNotFoundError('Compartment %s not found in current list' %(compartmentIds[cindex]))
                    
                    # Build the TemplateReaction.        
                    reaction['id'] = '%s_%s' %(reactionId, idcomp) # Use first compartment for suffix
                    reaction['name'] = masterReaction['name']
                    reaction['deltaG'] = masterReaction['deltag']
                    reaction['deltaGErr'] = masterReaction['deltagerr']
                    reaction['status'] = masterReaction['status']
                    reaction['reversibility'] = masterReaction['reversibility']
                    reaction['direction'] = fields[fieldNames['direction']]
                    if fields[fieldNames['gfdir']] == 'null':
                        reaction['GapfillDirection'] = ''
                    else:
                        reaction['GapfillDirection'] = fields[fieldNames['gfdir']]
                    reaction['type'] = fields[fieldNames['type']]
                    reaction['maxforflux'] = float(100)
                    reaction['maxrevflux'] = float(-100)
                    reaction['templatecompartment_ref'] = '~/compartments/id/'+compartmentIds[0]
                    reaction['base_cost'] = float(fields[fieldNames['base_cost']])
                    reaction['forward_penalty'] = float(fields[fieldNames['forward_cost']])
                    reaction['reverse_penalty'] = float(fields[fieldNames['reverse_cost']])
                    reaction['templateReactionReagents'] = list()
                    # Stoichiometry format is n:cpdid:c:i:"cpdname"
                    if len(masterReaction['stoichiometry']) > 0:
                        try:
                            reagents = masterReaction['stoichiometry'].split(';')
                            for rindex in range(len(reagents)):
                                parts = reagents[rindex].split(':')
                                compartmentIndex = int(parts[2])
                                compCompound = self.addCompCompound(parts[1], compartmentIds[compartmentIndex])
                                templateReactionReagent = dict()
                                templateReactionReagent['templatecompcompound_ref'] = '~/compcompounds/id/'+compCompound['id']
                                templateReactionReagent['coefficient'] = float(parts[0])
                                reaction['templateReactionReagents'].append(templateReactionReagent)
                        except IndexError as e:
                            raise ReactionFormatError('Reaction %s on line %d has invalid stoichiometry "%s" in master reaction' %(reaction['id'], linenum, masterReaction['stoichiometry']))
                    reaction['templatecomplex_refs'] = list()
                    if reaction['type'] == 'conditional' and fields[fieldNames['complexes']] == 'null':
                        raise NoComplexesError('Reaction %s is of type conditional and no complexes are specified' %(reactionId))
                    if fields[fieldNames['complexes']] != 'null':
                        complexes = fields[fieldNames['complexes']].split('|')
                        for cindex in range(len(complexes)):
                            if complexes[cindex] in self.complexes:
                                reaction['templatecomplex_refs'].append('~/complexes/id/'+complexes[cindex])
                            else:
#                                print 'Reaction %s on line %d refers to complex %s which is not found' %(reaction['id'], linenum, complexes[cindex])
                                raise ComplexNotFoundError('Reaction %s on line %d refers to complex %s which is not found' %(reaction['id'], linenum, complexes[cindex]))
                        
                if includeLinenum:
                    reaction['linenum'] = linenum

                # Check for duplicates.
                if reaction['id'] not in self.reactions:
                    self.reactions[reaction['id']] = reaction
                else:
                    raise DuplicateReactionError('Reaction %s on line %d is a duplicate' %(reaction['id'], linenum))
        
        return

    def addCompCompound(self, compoundId, compartmentId):
        ''' Add a compound to the model template.
        
            @param compoundId: Compound ID in master compound list
            @param compartmentId: Compartment ID
            @return TemplateCompCompound structure
        '''
        
        # Get the compound from the master list.
        try:
            masterCompound = self.masterCompoundsList[self.masterCompounds[compoundId]]
        except KeyError as e:
            raise CompoundNotFoundError('Compound %s not found in master biochemistry' %(compoundId))

        # Check for obsolete compound.
        if masterCompound['is_obsolete']:
            if masterCompound['linked_compound'] != 'null':
                linkIds = masterCompound['linked_compound'].split(';')
                # One of the compounds in the list is not obsolete.
                for index in range(len(linkIds)):
                    try:
                        linkCompound = self.masterCompoundsList[self.masterCompounds[linkIds[index]]]
                        if not linkCompound['is_obsolete']:
                            print 'NOTICE: Obsolete compound %s replaced by %s' %(masterCompound['id'], linkCompound['id'])
                            masterCompound = linkCompound
                            break
                    except KeyError as e:
                        raise ObsoleteCompoundError('Compound %s is obsolete and replacement %s not found' %(masterCompound['id'], linkIds[index]))
                if masterCompound['is_obsolete']:
                    raise ObsoleteCompoundError('Compound %s is obsolete and all replacements are obsolete' %(masterCompound['id']))
            else:
                raise ObsoleteCompoundError('Compound %s is obsolete and no replacement is specified' %(masterCompound['id']))
        
        # If needed, create a new TemplateCompound and add it to the Model Template.
        if compoundId not in self.compounds:
            try:
                compound = dict()
                compound['id'] = masterCompound['id']
                compound['name'] = masterCompound['name']
                compound['abbreviation'] = masterCompound['abbreviation']
                compound['isCofactor'] = masterCompound['is_cofactor']
                compound['aliases'] = masterCompound['aliases']
                compound['defaultCharge'] = masterCompound['charge']
                compound['mass'] = masterCompound['mass']
                if compound['mass'] == 'null':
                    compound['mass'] = 0
                compound['deltaG'] = masterCompound['deltag']
                compound['deltaGErr'] = masterCompound['deltagerr']
                compound['formula'] = masterCompound['formula']
                self.compounds[compound['id']] = compound
            except KeyError as e:
                raise CompoundKeyError('Missing key in compound %s: %s' %(masterCompound['id'], e.message))
        else:
            compound = self.compounds[compoundId]

        # The id for TemplateCompCompound includes the compartment suffix.
        try:
            compartment = self.compartments[compartmentId]
        except KeyError as e:
            raise CompartmentNotFoundError('Compartment %s not found in current list' %(compartmentId))
        id = '%s_%s' %(compoundId, compartmentId)
        
        # If needed, create a new TemplateCompCompound and add it to the Model Template.
        if id not in self.compCompounds:
            compCompound = dict()
            compCompound['id'] = id
            compCompound['templatecompound_ref'] = '~/compounds/id/'+compoundId
            compCompound['charge'] = compound['defaultCharge'] # @todo Not sure how charge could be different
            if compartment['id'] == 'e':
                compCompound['maxuptake'] = 100.0 # Set a maximum for the extracellular compartment
            else:
                compCompound['maxuptake'] = 0.0
            compCompound['templatecompartment_ref'] = '~/compartments/id/'+compartmentId
            self.compCompounds[id] = compCompound
        else:
            compCompound = self.compCompounds[id]

        return compCompound
