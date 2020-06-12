# ModelSEED Biochemistry Compounds

The compounds are found in the compounds.tsv and compounds.json
files. They are described here:

1. **id**: Unique ID for the compound in the format cpdNNNNN where NNNNN is a five digit number (e.g. cpd00001)
2. **abbreviation**: Short name of compound
3. **name**: Long name of compound
4. **formula**: Standard chemical format (using Hill system) in protonated form to match reported charge
5. **mass**: Mass of compound or "null" when unknown
6. **source**: Source database of compound (currently only source is ModelSEED)
7. **inchikey**: Structure of compound using IUPAC International Chemical Identifier (InChI) format
8. **charge**: Electric charge of compound
9. **is_core**: True if compound is in core biochemistry (currently all compounds are set to true)
10. **is_obsolete**: True if compound is obsolete and replaced by different compound (currently all compounds are set to false)
11. **linked_compound**: List of compound IDs separated by semicolon related to this compound or "null" if not specified (used to link an obsolete compound to replacement compound)
12. **is_cofactor**: True if compound is a cofactor (currently all compounds are set to false)
13. **deltag**: Value for change in free energy of compound or "null" when unknown
14. **deltagerr**: Value for change in free energy error of compound or "null" when unknown
15. **pka**: Acid dissociation constants of compound (see below for description of format)
16. **pkb**: Base dissociation constants of compound (see below for description of format)
17. **abstract_compound**: True if compound is an abstraction of a chemical concept (currently all compounds are set to null)
18. **comprised_of**:  or "null" if not specified (currently all compounds are set to null)
19. **aliases**: List of alternative names of compound separated by semicolon or "null" if not specified (see below for description of format)
20. **smiles**: Structure of compound using Simplified Molecular-Input Line-Entry System (SMILES) format
21. **notes**: Abbreviated notation used to store derived information about the compound.

### Format of pka and pkb
The pka and pkb fields are in this format:

    fragment:atom:value

"fragment" is the index of the molecular fragment. This is almost
always 1 but there are a few molecular structures that contain more
than one distinct structure. "atom" is the index of the atom and
"value" is the dissociation constant value.  The pkas or pkbs of
multiple atoms are separated by a semicolon.  For example, this is the
pka for NAD:

    1:17:1.8;1:18:2.56;1:6:12.32;1:25:11.56;1:35:13.12

### Format of aliases
An alias is in this format:

   "source:value"
	
where "source" is the name of the alternative database and "value" is
the name or ID in the alternative database. Multiple aliases are
separated by a semicolon.  For example, this is the list of aliases
for Cobamide (cpd00181):

    "KEGG:C00210";"name:Cobamide";"searchname:cobamide";"ModelSEED:cpd00181";"KBase:kb|cpd.181"

### Notes abbreviations

### Sources