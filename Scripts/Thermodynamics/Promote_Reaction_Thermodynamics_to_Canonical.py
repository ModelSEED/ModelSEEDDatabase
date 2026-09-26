#!/usr/bin/env python
"""RETIRED 2026-09-11. Does nothing; kept so any caller fails loudly.

This promoted the best per-source estimate into the canonical top-level
`deltag` / `deltagerr` fields. Those fields were removed: energies now ship
per source under `thermodynamics`, and the single recommended direction ships
in `reversibility`, written by
Scripts/Thermodynamics/Apply_Evidence_Grades_And_Recommendation.py.

There is no longer a canonical energy to promote to.
"""
import sys

print(__doc__)
sys.exit(0)
