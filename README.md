[![Build Status](https://travis-ci.org/ModelSEED/ModelSEEDDatabase.svg?branch=master)](https://travis-ci.org/ModelSEED/ModelSEEDDatabase)
# ModelSEEDDatabase
This repository contains the definitive copy of the biochemistry and metadata used to construct models using the ModelSEED/ProbAnno approach

You can contribute by submitting revisions to the biochemistry though pull requests

# Python Environment

1) Install Conda:
https://conda.io/projects/conda/en/latest/user-guide/install/index.html
2) Set up RDKit environment with python 3
conda create -c rdkit -n rdkit-env rdkit python=3
3) Activate:
conda activate rdkit-env
4) Install openbabel:
conda install -c openbabel openbabel
5) Export path to local python libraries
export PYTHONPATH=/Users/seaver/Projects/ModelSEEDDatabase/Libs/Python/BiochemPy/
6) Test it by running:
./Scripts/Structures/Print_Structure_Formula_Charge.py
and
git status -s
(there should be no change in the formula and charge in the compounds.tsv file, if there is, our recommendation is that you panic and throw your laptop out of the window.)