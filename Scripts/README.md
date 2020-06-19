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