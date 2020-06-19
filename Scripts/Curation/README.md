# ModelSEED Biochemistry Database Scripts (Curation)

We developed two set of scripts to help us consistently deal with two kinds of problems that we encountered in building the database:

* [Redundancy](Redundancies)

  In this category, we simply find that there are two separate
  compounds that actually represent the same identical molecular
  structure and, for whatever reason, were not merged when building
  the database. Here we merge them. Doing so involves merging all the
  attached data for both of them, and rendering one of them
  obsolete for the purposes of metabolic modeling.

* [Conflict](Disambiguation/README.md)

  In this category, we find that two compounds from the external
  sources that were used to build the database were incorrectly
  merged, for whatever reason (a common reason is a generic
  synonym). To curate this problem, we wrote scripts that separated
  out the two original compounds into two separate entities, and
  sorted out the associated data, espcially aliases and
  structures. Sometimes, this means creating an entirely new compound,
  and sometimes, the second compound already exists, in which case
  we'd merge the associated data with it. The process of separating
  the data out, using provenance where we can is convoluted and we
  describe it a little more in the Disambiguation directory.

You feed a single compound identifier into 

Print_Compound_to_Disambiguate.py

eg:

./Print_Compound_to_Disambiguate.py cpd03657

You edit the resulting file in Objects, i.e.

Objects/cpd03657_Sam_Seaver_Object.json

(it uses your git name)

To edit, you have to change any value from true to false if it doesn't belong to the originating compound
(you might have to create a new compound, though the code attempts to find one in the current biochemistry)

Then you re-integrate the objet in

./Integrate_Disambiguated_Compound.py Objects/cpd03657_Sam_Seaver_Object.json

If it finds that the disambiguation affects any reactions, then you get a reactions object
which you have to edit, then you feed it back in:

./Integrate_Disambiguated_Reactions.py <example>

Finally, you have to re-run code for updating structures, formulas, aliases, reaction balancing, and proton balancing

cd ../Structures/
./List_ModelSEED_Structures.py
./Update_Compound_Structures_Formulas_Charge.py

cd ../Biochemistry/
./Update_Compound_Aliases.py
./Rebalance_Reactions.py
./Adjust_Reaction_Protons.py

Use git diff to check everything, and commit and make a new PR.
Be sure to add the disambiguation objects in the Objects directory too