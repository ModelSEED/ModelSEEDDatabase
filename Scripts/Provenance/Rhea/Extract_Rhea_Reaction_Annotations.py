#!/usr/bin/env python
import rdflib
from rdflib import Graph

graph = Graph()
graph.parse('Data/rhea.rdf')

url_list = ('purl.uniprot.org', 'identifiers.org', 'purl.obolibrary.org', 'rdf.ncbi.nlm.nih.gov')

for subject, predicate, object in graph.triples((None, None, None)):

	is_found=False
	for url in url_list:
          if(url in object and 'Compound' not in subject):
               #print(str(subject),str(predicate),str(object))
               rhea = str(subject).split('/')[-1]
               (db,id) = str(object).split('/')[-2:]
               print('\t'.join([rhea,db,id]))
               #is_found=True
          
	if is_found is True:
          break
