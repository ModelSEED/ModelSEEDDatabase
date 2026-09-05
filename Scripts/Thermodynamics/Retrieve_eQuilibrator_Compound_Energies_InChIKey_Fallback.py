#!/usr/bin/env python
"""Structure-based InChIKey fallback for eQuilibrator compound-energy retrieval.

The default Retrieve_eQuilibrator_Compound_Energies.py uses MetaNetX IDs
only. That flow misses two classes of ModelSEED compound: those with no
MetaNetX ID at all, and those whose MetaNetX ID isn't in eQuilibrator's
cache. This script fills both by asking eQuilibrator directly with the
compound's InChIKey.

Per Sam's constraint (2026-08-04): fallbacks MUST be structure-based
(InChIKey or canonical SMILES) to avoid identifier-remap risk. InChIKey
is a hash of the InChI, so it satisfies that constraint.

Writes:
  Biochemistry/Thermodynamics/eQuilibrator/
    InChIKey_Fallback_Compound_Energies.tbl

Columns: inchikey, dG0_prime_kcal, uncertainty_kcal, source_cpd
"""

if __name__ == "__main__":
    # Validate arguments BEFORE importing anything or touching the database.
    # These scripts mutate the database, and without this an unknown flag or a
    # mistyped mode was silently ignored and the script ran with its defaults:
    # asking Estimate_Reaction_Reversibility.py for --help rewrote 122 files.
    # Placed above the imports so --help works even where a dependency is
    # missing from the path.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import os
import sys
import json
import glob
from equilibrator_api import ComponentContribution, Q_, Reaction

# --- Setup ---
cc = ComponentContribution()
cc.p_h = Q_(7.0)
cc.ionic_strength = Q_("0.25M")
cc.temperature = Q_("298.15K")

REPO_ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."
COMPOUND_JSON_GLOB = os.path.join(REPO_ROOT, 'Biochemistry', 'compound_*.json')
OUTPUT = os.path.join(REPO_ROOT, 'Biochemistry', 'Thermodynamics',
                      'eQuilibrator',
                      'InChIKey_Fallback_Compound_Energies.tbl')

# --- Collect candidate compounds ---
candidates = []
for path in sorted(glob.glob(COMPOUND_JSON_GLOB)):
    with open(path) as fh:
        for c in json.load(fh):
            formula = c.get('formula') or ''
            smiles = c.get('smiles') or ''
            inchikey = c.get('inchikey') or ''
            td = c.get('thermodynamics')
            has_eq = isinstance(td, dict) and 'eQuilibrator' in td
            if smiles and 'R' not in formula and inchikey and not has_eq:
                candidates.append((c['id'], inchikey))

sys.stderr.write(
    f"Attempting InChIKey lookup + ΔfG computation for "
    f"{len(candidates):,} ModelSEED compounds without an existing "
    f"'eQuilibrator' entry.\n"
)

# --- Query eQuilibrator ---
# Correct API: cc.get_compound(inchi_key) returns the Compound directly,
# then a formation reaction Reaction({compound: 1}) lets us call
# cc.standard_dg_prime() to get ΔfG'.
os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
found_exact = 0
found_connectivity = 0
not_in_cache = 0
in_cache_no_dg = 0
with open(OUTPUT, 'w') as out_fh:
    out_fh.write("inchikey\tdG0_prime_kcal\tuncertainty_kcal\tsource_cpd\tmatch_type\n")
    for i, (cpd_id, inchikey) in enumerate(candidates):
        if i and i % 500 == 0:
            sys.stderr.write(
                f"  progress: {i:,} / {len(candidates):,}  "
                f"(exact: {found_exact:,}, connectivity: {found_connectivity:,}, "
                f"not in cache: {not_in_cache:,})\n"
            )
        # Priority 1: exact InChIKey match — same InChI (skeleton + stereo +
        # protonation). Zero ambiguity.
        matches = cc.ccache.search_compound_by_inchi_key(inchikey)
        match_type = 'exact'
        # Priority 2: connectivity-block match — same skeleton, different
        # stereo/protonation. When multiple, pick the lowest-id (canonical/
        # oldest registration). Tagged as 'connectivity' so downstream can
        # decide whether to use.
        if not matches:
            block = inchikey.split('-', 1)[0]
            matches = cc.ccache.search_compound_by_inchi_key(block)
            match_type = 'connectivity'
        if not matches:
            not_in_cache += 1
            continue
        compound = min(matches, key=lambda c: c.id)
        try:
            rxn = Reaction({compound: 1})
            result = cc.standard_dg_prime(rxn)
            dg = str(result.value.to('kilocal / mole').magnitude)
            unc = str(result.error.to('kilocal / mole').magnitude)
            out_fh.write('\t'.join([inchikey, dg, unc, cpd_id, match_type]) + '\n')
            if match_type == 'exact':
                found_exact += 1
            else:
                found_connectivity += 1
        except Exception:
            in_cache_no_dg += 1
            continue

sys.stderr.write(
    f"\nDone.\n"
    f"  Exact InChIKey match + ΔfG:        {found_exact:,}\n"
    f"  Connectivity-block match + ΔfG:    {found_connectivity:,}\n"
    f"    Total recovered:                 {found_exact + found_connectivity:,}\n"
    f"  Not in cache under any match:      {not_in_cache:,}\n"
    f"  In cache but ΔfG failed:           {in_cache_no_dg:,}\n"
    f"Output: {OUTPUT}\n"
)
