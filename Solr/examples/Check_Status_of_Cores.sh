#!/bin/bash

# checking status of loaded cores
curl http://localhost:8983/solr/admin/cores?action=STATUS
