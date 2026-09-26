For each of the primary sources that we use for integrating with the database, we don't keep their files
within this repository, but we have scripts that we use to extract the necessary information from
their original files. We do keep the derived and re-formatted data in ../../Biochemistry/Provenance
along with the report of how each compound and reaction was merged with the database.

The scripts for processing the data from each source is in each folder, and we have some database-specific notes:

For Rhea:

We downloaded the entire database as a single file in RDF format.

They use three sources/types of compounds, the ones that are registered with ChEBI, using ChEBI identifiers, 
these, for the most part, have a distinct structure and an accompanying mol file, GENERIC types that have one 
or more R groups and apparently are not found in ChEBI, but are using in metabolic modeling as abstract 
compounds to "complete" reactions or pathways that aren't fully understood, these ones don't have molecular 
structures associated with them, and finally POLYMER types with SRUs of an indeterminate length. These have
molecular structures but because of the way SRUs are encoded in the mol format, the indeterimnate size is not
encoded in the resulting InChI so we don't use them (we match names only with the `-no` flag)

As we extracted ChEBI identifiers for most of the Rhea compounds from the RDF file, but Rhea also uses their 
own identifiers for in-house
(from: Scripts/Biochemistry/Update)
./Add_New_Compounds.py ../../../Biochemistry/Provenance/Rhea/ChEBI_ID_Name_InChIKey.tsv ChEBI -r -s
./Add_New_Compounds.py ../../../Biochemistry/Provenance/Rhea/Rhea_Generic_Polymer_Names.tsv Rhea -r -s -no
./Add_New_Reactions.py ../../../Biochemistry/Provenance/Rhea/Rhea_reactions.tsv ChEBI,Rhea Rhea -r -s

For KEGG:

We downloaded the original compound and reaction files in KEGG format
 
./Add_New_Compounds.py ../../../Biochemistry/Provenance/KEGG/KEGG_compounds.tsv KEGG -r -s
./Add_New_Reactions.py ../../../Biochemistry/Provenance/KEGG/KEGG_reactions.tsv KEGG KEGG -r -s

For MetaCyc:

We downloaded the compounds.dat and reactions.dat files in PathwayTools format

./Add_New_Compounds.py ../../../Biochemistry/Provenance/MetaCyc/MetaCyc_compounds.tsv MetaCyc -r -s
./Add_New_Reactions.py ../../../Biochemistry/Provenance/MetaCyc/MetaCyc_reactions.tsv MetaCyc MetaCyc -r -s
