# ModelSEED Biochemistry Database Scripts (Curation)

We developed two set of scripts to help us consistently deal with two kinds of problems that we encountered in building the database:

# [Redundancy](Redundancies)

  In this category, we simply find that there are two separate
  compounds that actually represent the same identical molecular
  structure and, for whatever reason, were not merged when building
  the database. Here we merge them. Doing so involves merging all the
  attached data for both of them, and rendering one of them
  obsolete for the purposes of metabolic modeling.

# [Conflict](Disambiguation/README.md)

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
