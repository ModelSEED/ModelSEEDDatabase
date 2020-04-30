#!/usr/bin/env python
from urllib.request import urlopen
import json
connection = urlopen('http://localhost:8983/solr/admin/cores?action=STATUS')
response = json.load(connection)
print("Number of compounds :",response['status']['compounds']['index']['numDocs']
print("Number of reactions :",response['status']['reactions']['index']['numDocs']

# Print the name of each document.

#for document in response['response']['docs']:
#  print("  Name =", document['name'])
