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

