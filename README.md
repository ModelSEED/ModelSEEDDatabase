[![Build Status](https://travis-ci.org/ModelSEED/ModelSEEDDatabase.svg?branch=master)](https://travis-ci.org/ModelSEED/ModelSEEDDatabase)
# ModelSEEDDatabase
This repository contains the definitive copy of the biochemistry and metadata used to construct models using the ModelSEED/ProbAnno approach

You can contribute by submitting revisions to the biochemistry though pull requests

# Python Environment

1) Install Conda:
https://conda.io/projects/conda/en/latest/user-guide/install/index.html

2) Set up a python 3 conda environment:
```
 conda create -n msd-env python=3
```
 3) Activate:
```
 conda activate rdkit-env
```

4) Install rdkit and openbabel for handling biochemical structures:
```
 conda install -c rdkit rdkit
 conda install -c openbabel openbabel
```

5) Install eQuilibrator for retrieving thermodynamics data:
```
conda install wxPython
pip install quilt
```

6) Export path to local python libraries
```
export PYTHONPATH=$PYTHONPATH:<path-to-repository>/ModelSEEDDatabase/Libs/Python/
```

7) Test it by running:
```
./Scripts/Structures/Print_Structure_Formula_Charge.py
git status -s
```
(there should be no change in the formula and charge in the compounds.tsv file, if there is, our recommendation is that you panic and throw your laptop out of the window.)