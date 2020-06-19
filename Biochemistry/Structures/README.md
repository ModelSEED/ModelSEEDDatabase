# ModelSEED Biochemistry Database Structures

Here we keep all the molecular structures we use for the ModelSEED
Biochemistry Database. We download and use molecular structures from
KEGG, MetaCyc and MetaNetX, each stored as SMILES, InChI, and InChIKey
representations in their own folder, and within the
`All_ModelSEED_Structures.txt` file.

We use Marvin to protonate each mol file at a pH of 7 cannot handle
every single file for various reasons, so all files that contain the
word "Original" were generated from the set of uncharged molecular
structure files that were downloaded from their respective sources,
and all files with the word "Charged" concern the files that were
sucessfully protonated by Marvin.

Marvin was also used to generate the pKas and pKbs of the atoms in
each molecular structure.

We also parsed the KEGG molecular structure files for SRUs (structural
repeating units). The information stored for these are passed over by
Marvin when protonating, and not repeated in the resulting output
structure. This caused problems when merging compounds that wre
clearly not the same compound, so we keep a record of which of the
KEGG molecular structure files contain these SRUs.