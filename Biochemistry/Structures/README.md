# ModelSEED Biochemistry Database Structures

Here we keep all the molecular structures we use for the ModelSEED Biochemistry Database


MarvinBeans:

	The version of MarvinBeans used to protonate each mol file at
	a pH of 7 cannot handle every single file for various reasons,
	so all files that contain the word "Original" were generated
	from the set of uncharged MOL files that were downloaded from
	their respective sources, and all files with the word
	"Charged" concern the files that were sucssfully protonated by
	MarvinBeans.

MolAnalysis: 

	ModelSEEDCore code was used to derive the chemical groups that
	compose the entire structure of a molecular and are in turn
	used for thermodynamic heuristics.

InChI:

	InChI strings cannot be generated from files containing
	"R" groups or pseudo-atoms, amongst other reasons, and as
	such, the number of structures represented in the
	"Original_InChI" files is less than the number of structures
	represented in the "Original_MolAnalysis" files.
