MarvinBeans:

	The version of MarvinBeans used to protonate each mol file at
	a pH of 7 cannot handle every single file for various reasons,
	so all files that contain the word "Original" were generated
	from the set of uncharged MOL files that were downloaded from
	their respective sources, and all files with the word
	"Charged" concern the files that were sucssfully protonated by
	MarvinBeans.

MolAnalysis: 

	ModelSEEDCore code was used to generate the set of groups
	that compose the chemical structure and are in turn used for
	thermodynamic heuristics.

InChI:

	InChI strings cannot be generated from files containing
	"R" groups or pseudo-atoms, amongst other reasons, and as
	such, the number of structures represented in the
	"Original_InChI" files is less than the number of structures
	represented in the "Original_MolAnalysis" files.
	
Search:

	Radicals and unusual valences not captured in the otherwise
	canonical InChI strings, but are printed out in the auxilliary
	information. In order to capture this, so that radicals will
	not be merged with their elementally-identical counter-parts,
	this information is appended to the InChI strings in the files
	containing the word "Search". It is the "Search" strings that
	I use to merge compound structures.

Counts: 

   16693 KEGG_Original_MolAnalysis.tbl
   16692 KEGG_Charged_MolAnalysis.tbl

   15450 KEGG_Original_InChI.txt
   15448 KEGG_Charged_InChI.txt
   15449 KEGG_Search_InChI.txt

   11181 MetaCyc_Original_MolAnalysis.tbl
   13698 MetaCyc_Charged_MolAnalysis.tbl

   11660 MetaCyc_Original_InChI.txt
   11659 MetaCyc_Charged_InChI.txt
   11659 MetaCyc_Search_InChI.txt