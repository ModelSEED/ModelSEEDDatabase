# ModelSEED Biochemistry Reactions

The reactions are found in the reactions.tsv and reactions.json
files. They are described here:

1. **id**: Unique ID for the reaction in the format rxnNNNNN where NNNNN is a five digit number (e.g. rxn03789)
2. **abbreviation**: Short name of reaction
3. **name**: Long name of reaction
4. **code**: Definition of reaction expressed using compound IDs and before protonation (see below for description of format)
5. **stoichiometry**: Definition of reaction expressed in stoichiometry format (see below for description of format)
6. **is_transport**: True if reaction is a transport reaction 
7. **equation**: Definition of reaction expressed using compound IDs and after protonation (see below for description of format)
8. **definition**: Definition of reaction expressed using compound names (see below for description of format)
9. **reversibility**: Reversibility of reaction where ">" means right directional, "<" means left directional, "=" means bi-directional, and "?" means unknown (_Need a better description of reversibility and direction_)
10. **direction**: Direction of reaction where ">" means right directional, "<" means left directional, and "=" means bi-directional
11. **abstract_reaction**: _Need definition_ or "null" if not specified (currently all reactions are set to null)
12. **pathways**: Pathways reaction is a part of or "null" if not specified (currently all reactions are set to null)
13. **aliases**: List of alternative names of reaction separated by semicolon or "null" if not specified (format is the same as Compounds file)
14. **ec_numbers**: Enzyme Commission numbers of enzymes that catalyze reaction or "null" if not specified (currently all reactions are set to null)
15. **deltag**: Value for change in free energy of reaction or 10000000 when unknown
16. **deltagerr**: Value for change in free energy error of reaction or 10000000 when unknown
17. **compound_ids**: List of compound IDs separated by semicolon for compounds involved in reaction
18. **status**: String describing status of the reaction with multiple values delimited with a "|" character.  See below for details.
19. **is_obsolete**: True if reaction is obsolete and replaced by different reaction
20. **linked_reaction**: List of reaction IDs separated by semicolon related to this reaction or "null" if not specified (used to link an obsolete reaction to replacement reaction)
21. **notes**: Abbreviated notation used to store derived information about the reaction
22. **source**: Source database of reaction (currently only source is ModelSEED)

### Format of stoichiometry field

### Format of reaction definition using compound IDs 
Each compound participating in the reaction is in this format:

	(n) cpdid[m]

where "n" is the compound coefficient, "cpdid" is the compound ID, and "m" is a compartment index number.  Compounds are separated by "+" and reactant and product compounds are delimited by the direction symbol.  For example, this is the definition of reaction rxn00001:

	(1) cpd00001[0] + (1) cpd00012[0] <=> (2) cpd00009[0] + (1) cpd00067[0]

### Format of reaction definition using compound names 
Each compound participating in the reaction is in this format:

	(n) cpdname[m]

where "n" is the compound coefficient, "cpdname" is the compound name, and "m" is a compartment index number.  Compounds are separated by "+" and reactant and product compounds are delimited by the direction symbol.  For example, this is the definition of reaction rxn00001:

	(1) H2O[0] + (1) PPi[0] <=>  (2) Phosphate[0] + (1) H+[0]

### Format of reaction stoichiometry
Each compound participating in the reaction is in this format:

	n:cpdid:m:i:"cpdname"

where "n" is the compound coefficient and a negative number indicates
a reactant and a positive number indicates a product, "cpdid" is the
compound ID, "m" is the compartment index number, "i" is the community
index number (I think this is no longer needed), and "cpdname" is the
compound name.  Compounds are separated by semicolon.  For example,
this is the stoichiometry of reaction rxn00001:

	-1:cpd00001:0:0:"H2O";-1:cpd00012:0:0:"PPi";2:cpd00009:0:0:"Phosphate";1:cpd00067:0:0:"H+"

### Reaction status values

The reaction status field is a string with one or more
values. Multiple values are separated with a "|" character.  There are
multiple values when a reaction is updated after mass and charge
balancing.  If "OK" is one of the values, the reaction is valid and
additional values describe the changes made to the reaction to balance
it.  The status values are:

* OK means the reaction is valid.  If "OK" is the only value, then the reaction was valid with no changes.
* MI means there is a mass imbalance. The remainder of the string after the first colon indicates what atoms are unbalanced and the number of atoms needed to balance the reaction.  Multiple atoms are separated by a "/" character.  A positive number means the righthand side of the reaction has more atoms and a negative number means the lefthand side of the reaction has more atoms.
* CI means there is a charge imbalance.  A positive number after the first colon means the righthand side of the reaction has a larger charge and a negative number means the lefthand side of the reaction has a larger charge.
* HB means the reaction has been balanced by adding hydrogen to it.
* EMPTY means reactants cancel out completely.
* CPDFORMERROR means at least one compound either has no formula or has an invalid formula.

For example, rxn00277 has this definition:

	(1) Glycine[0] <=> (1) HCN[0] or
	(1) C2H5NO2 <=> (1) CHN
	
and its status shows that it has a mass imbalance with 1 extra carbon, 4 extra hydrogen, 2 extra oxygen atoms on the lefthand side of the reaction:

	MI:C:-1/H:-4/O:-2

And rxn00008 has this definition:

	(2) H2O[0] <=> (1) H2O2[0] + (2) H+[0]

and its status shows it has a charge imbalance with the righthand side of the reaction having a larger charge:

	CI:2
