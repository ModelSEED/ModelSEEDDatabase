#!/usr/bin/env python
from urllib.request import urlopen
import json

SOLR_URL='https://modelseed.org'

## Uncomment the line below
## if you would like to query your local Solr server
#SOLR_URL='http://localhost:8983'

# Retrieving compounds with 'water' in their synonyms
connection = urlopen(SOLR_URL+'/solr/compounds/select?wt=json&q=aliases:water&fl=name,id,formula,charge,aliases')
response = json.load(connection)
for document in response['response']['docs']:
  print(document['name'],document['id'],document['formula'],document['charge'],document['aliases'])

# Retrieving reactions with 'water' in their synonyms
connection = urlopen(SOLR_URL+'/solr/reactions/select?wt=json&q=aliases:water&fl=name,id,definition')
response = json.load(connection)
for document in response['response']['docs']:
  print(document['name'],document['id'],document['definition'])
