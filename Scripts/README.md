# ModelSEED Biochemistry Database Script Environment

We've developed a small python library for loading and manipulating
the biochemistry database according to our needs and we written a set
of scripts that uses the library for our needs.

The main set of scripts that we used for the most recent publication of the ModelSEED Biochemistry Database are:

* [Biochemistry](Biochemistry/README.md)

* [Structures](Structures/README.md)

* [Thermodynamics](Thermodynamics/README.md)

* [Curation](Curation/README.md)

* [Statistics](Statistics/README.md)

In the last case, the scripts were used to generate some of the figures and tables in our latest publication.

### Python Environment

All these scripts run in a Python 3 environment, and some require additional packages which we explain here:

1) Install Conda
https://conda.io/projects/conda/en/latest/user-guide/install/index.html

2) Set up a python 3 conda environment
```
 conda create -n msd-env python=3
```
 3) Activate
```
 conda activate msd-env
```

4) Install rdkit and openbabel for handling biochemical structures
```
 conda install -c rdkit rdkit
 conda install -c openbabel openbabel
```

5) Install eQuilibrator for retrieving thermodynamics data
```
conda install wxPython
pip install quilt
pip install equilibrator_api
```
_(you might not need to install `wxPython`, I had trouble with dependencies on a mac)_

6) Export path to local python libraries
```
export PYTHONPATH=$PYTHONPATH:<path-to-repository>/ModelSEEDDatabase/Libs/Python/
```

7) Test
```
./Biochemistry/Reprint_Biochemistry.py
git status -s
```
(The script should run without throwing any errors and there should be
no change in the biochemistry data)

## The argument guard

Most scripts here open with this, above the imports:

```python
if __name__ == "__main__":
    # Argument guard -- see "The argument guard" in Scripts/README.md.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()
```

**Why it exists.** These scripts mutate the database. Without the guard an
unknown flag or a mistyped mode was silently ignored and the script ran with its
defaults: asking `Estimate_Reaction_Reversibility.py` for `--help` rewrote 122
files. The guard turns that into an `argparse` error before anything is read or
written.

**Why it sits above the imports.** So `--help` still works on a machine where one
of the script's dependencies is missing from the path. Move it below the imports
and `--help` fails with an ImportError instead of printing the docstring.

**Keep it first.** Anything placed before it -- an import with a side effect, a
module-level path lookup, a database open -- runs before the arguments have been
checked, which is the failure the guard exists to prevent.

Scripts that take real arguments extend the same block with their own
`add_argument` calls; `Scripts/Structures/Run_Marvin_pKas.py` is an example.
