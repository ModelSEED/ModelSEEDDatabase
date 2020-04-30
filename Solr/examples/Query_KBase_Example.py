#!/usr/bin/env python
# pip install --extra-index-url https://pypi.anaconda.org/kbase/simple releng-client
from relation_engine_client import REClient
re_client = REClient("https://kbase.us/services/relation_engine_api")

# Retrieving compounds with 'water' in their synonyms
response = re_client.stored_query("search_compounds",{"search_text":"water"})
for document in response['results']:
  print(document['name'],document['id'],document['formula'],document['charge'],document['aliases'])

# Retrieving reactions with 'water' in their synonyms
response = re_client.stored_query("search_reactions",{"search_text":"water"})
for document in response['results']:
  print(document['name'],document['id'],document['definition'])
