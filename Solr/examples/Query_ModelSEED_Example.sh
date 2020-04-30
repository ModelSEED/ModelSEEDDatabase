#!/bin/bash
# Retrieving compounds with 'water' in their synonyms
curl "https://modelseed.org/solr/compounds/select?wt=csv&q=aliases:water&fl=name,id,formula,charge,aliases"

# Retrieving reactions with 'water' in their synonyms
curl "https://modelseed.org/solr/reactions/select?wt=csv&q=aliases:water&fl=name,id,definition"

