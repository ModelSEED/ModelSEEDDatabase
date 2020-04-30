#!/bin/bash
URL="https://kbase.us/services/relation_engine_api/api/v1/"
# Retrieving compounds with 'water' in their synonyms
curl -X POST "${URL}query_results?stored_query=search_compounds" -d '{ "search_text": "water" }'

# Retrieving reactions with 'water' in their synonyms
curl -X POST "${URL}query_results?stored_query=search_reactions" -d '{ "search_text": "water" }'
