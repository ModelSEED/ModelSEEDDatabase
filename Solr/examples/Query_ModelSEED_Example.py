#!/usr/bin/env python
from urllib.request import urlopen
import json

# Retrieving compounds with 'water' in their synonyms
connection = urlopen('https://modelseed.org/solr/compounds/select?wt=json&q=aliases:water&fl=name,id,formula,charge,aliases')
response = json.load(connection)
for document in response['response']['docs']:
  print(document['name'],document['id'],document['formula'],document['charge'],document['aliases'])

# Retrieving reactions with 'water' in their synonyms
connection = urlopen('https://modelseed.org/solr/reactions/select?wt=json&q=aliases:water&fl=name,id,definition')
response = json.load(connection)
for document in response['response']['docs']:
  print(document['name'],document['id'],document['definition'])
