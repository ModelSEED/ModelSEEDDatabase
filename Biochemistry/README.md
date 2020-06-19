# ModelSEED Biochemistry Database

The ModelSEED Biochemistry Database is found in this directory.  

There are quite a few files kept together in this directory and its
sub-directories, and we will describe here the most relevant ones.

The main biochemical data that we use, visualize, and serve, are in
this directory and consists of two sets of files, for compounds and
reactions, in two different formats:

* `compounds.tsv`
* `compounds.json`
* `reactions.tsv`
* `reactions.json`

We cover the fields we use to describe the biochemistry in the
[COMPOUNDS.md](COMPOUNDS.md) and [REACTIONS.md](REACTIONS.md) files respectively.

### Do Not Edit These Files Directly

The TSV files are in currently the _master_ files, and the JSON files
are created from the data in the TSV files on the fly. This is an
important point because, at this present time, if you edit the JSON
files, and then use our scripts to manage the database, you will lose
your changes. We'd like to make the point that users should use our
scripts, or write custom scripts that uses our python library to
interface with the data in the TSV files to avoid any confusion.  

All the scripts that we use to manage our data can be found in the
`../Scripts` directory, and there are many examples of what you can do
with the data. The underlying python library, that loads the data as
dictionaries for manipulation, can be found in
`../Libs/Python/BiochemPy`.

## Additional Data

We have collected and grouped data that underlies our construction of
the ModelSEED Biochemistry Database into these folders

# [Aliases](Aliases/README.md)

All of the ModelSEED Biochemistry came from external sources, the
majority of which comes from <a href="https://www.kegg.jp/">KEGG</a>
and <a href="https://metacyc.org/">MetaCyc</a>. Here, in the
[Aliases](Aliases) directory, we link every compound and reaction to
their sources.

# [Structures](Structures/README.md)

We attempt to link and consolidate our compounds around biochemical
structures, and in turn match reactions based on these structures. We
also use biochemical structures for computing compound and reaction
energies. In this way, we attempt to ground the core of our
biochemistry on actual physical chemistry. The structures we use are
downloaded from KEGG and MetaCyc and their derivatives are stored in
this directory.

# [Thermodynamics](Thermodynamics)

The results from our two different approaches to computing the
compound and reaction energies are stored in this directory.

# [Pathways](Pathways)

We aim to keep a list of subsystems from the ModelSEED environment as
well as external sources in this directory, but we have not yet begun
to maintain it.