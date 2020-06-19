#!/bin/bash
SOLR_URL="https://modelseed.org"

## Uncomment the line below
## if you would like to query your local Solr server
#SOLR_URL="http://localhost:8983"

# Retrieving compounds with 'water' in their synonyms
curl "${SOLR_URL}/solr/compounds/select?wt=csv&q=aliases:water&fl=name,id,formula,charge,aliases"

# Retrieving reactions with 'water' in their synonyms
curl "${SOLR_URL}/solr/reactions/select?wt=csv&q=aliases:water&fl=name,id,definition"

