#!/usr/bin/env python
## To install the KBase python client run:
## pip install --extra-index-url https://pypi.anaconda.org/kbase/simple releng-client
## The specification (DB schema) for the compounds and reactions are found here:
## https://github.com/kbase/relation_engine_spec/blob/master/stored_queries/
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
