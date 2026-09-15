#!/usr/bin/env python
"""RETIRED 2026-09-11. Does nothing; kept so the cascade fails loudly.

This script existed to pick one Marvin value per compound from the per-database
cascade (KEGG > MetaCyc > ChEBI > Rhea) and write it into the flat `pka` / `pkb`
fields. Those fields were removed from the compound records: protonation now
ships per tool under `pkas`, written by Add_Compound_pKa_Sources.py, which runs
the same cascade for the Marvin entry and additionally carries MolGpKa and
Literature.

Nothing was lost. The flat fields held the cascade-winning Marvin value, which
is exactly `pkas["Marvin"]`.
"""
import sys

print(__doc__)
print("Nothing to do. Use Scripts/Structures/Add_Compound_pKa_Sources.py.")
sys.exit(0)
