''' Exception thrown when source file header is invalid '''
class BadHeaderError(Exception):
    pass

''' Helper methods for working with common data structures. '''

class BaseHelper:
    
    def __init__(self):
        ''' Initialize object.
        '''

        pass

    def buildDictFromListOfObjects(self, objectList, key='id'):
        ''' Build a dictionary with the specified key from a list of objects.
    
            The value of each element in the returned dictionary is the complete
            object.  This copies the objects into a separate data structure.
    
            @param objectList: List of objects
            @param key: Field in object to use for keys in returned dictionary
            @return Dictionary mapping key to object
        '''
    
        objectDict = dict()
        for index in range(len(objectList)):
            obj = objectList[index]
            if obj is not None:
                objectDict[obj[key]] = obj
        return objectDict

    def buildIndexDictFromListOfObjects(self, objectList, key='id'):
        ''' Build a dictionary with the specified key to the index in the list.
    
            The value of each element in the returned dictionary is the index into
            the list of the object.  This allows for lookup of objects by the key
            without copying the objects.
    
            @param objectList: List of objects
            @param key: Field in object to use for keys in returned dictionary
            @return Dictionary mapping key to index of element in objectList
        '''
    
        indexDict = dict()
        for index in range(len(objectList)):
            obj = objectList[index]
            if obj is not None:
                indexDict[obj[key]] = index
        return indexDict
    
    def validateHeader(self, nameList, required):
        ''' Validate the header line in a source file.
        
            @param nameList: List of field names from first line of source file
            @param required: Dictionary of required field names
            @return Dictionary keyed by field names where value is index of field in nameList
            @raise BadHeaderError: Required field missing from header
        '''

        fieldNames = dict()
        for index in range(len(nameList)):
            fieldNames[nameList[index]] = index
        for req in required:
            if req not in fieldNames:
                raise BadHeaderError('Required field %s is missing from header' %(req))
        return fieldNames

    def addToList(self, sourceString, delim, destList):
        ''' Add to a list of items from a source string.
        
            @param sourceString: Source string with items separated by a delimiter
            @param delim: Delimiter separating items
            @param destList: List that items are appended to
            @return Nothing
        '''

        items = sourceString.split(delim)
        for index in range(len(items)):
            destList.append(items[index])
        return
