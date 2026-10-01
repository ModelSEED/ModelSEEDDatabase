#!/usr/bin/env python
"""Write Group-Contribution reaction energies by summing the per-compound
energies stored under the same label.

A reaction gets a real energy only when *every* reagent has a usable
per-compound energy (the original 'all-or-nothing' GC reaction rule);
otherwise the default sentinel is written. Obsolete-compound linking is
applied while building the eligible-compound set."""
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds, Reactions
import _thermo_helpers as th

LABEL = 'Group contribution'

th.run_reaction_aggregation_update(
    Reactions(), Compounds(), LABEL)
