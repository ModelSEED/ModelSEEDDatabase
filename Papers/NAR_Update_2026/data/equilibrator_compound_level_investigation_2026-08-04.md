# eQuilibrator compound-level ΔfG — investigation notes

**Date:** 2026-08-04
**Context:** Follow-up to the eQuilibrator rerun (Andrew's `equilibrator-rerun-modelseed-structures` branch, merged into dev at c1149b76) and the compound-level propagation (fb2bef4f). Aim: identify why coverage gaps remain and what a fix would look like.

## The gap after the propagation

- Non-R structured ModelSEED compounds: **30,479**
- With any deltag (per-source or top-level): **27,586** (91%)
- **Remaining gap: 2,893** structured compounds still without ΔfG

Of the compounds already covered by Andrew's rerun, 21,840 have an explicit `eQuilibrator` entry in their per-source thermodynamics dict. The rest of the coverage comes from historical top-level `deltag` values (Group Contribution and earlier eQuilibrator runs).

## The `Retrieve_eQuilibrator_Compound_Energies.py` script uses MetaNetX IDs only

Andrew's script constructs a per-compound query as:

```python
Reaction.parse_formula(equilibrator_calculator.get_compound, ' = ' + mnx)
```

where `mnx` is a MetaNetX ID from `Biochemistry/Structures/MetaNetX/Structures_in_ModelSEED_and_eQuilibrator.txt`. That file is the **intersection** of ModelSEED and eQuilibrator by structure (matching by InChIKey via `Find_eQuilibrator_Structures_in_ModelSEED.py`) but only carries MetaNetX-ID accessions. So the retrieval flow is:

1. Structure-based filter to find compounds present in both DBs (InChIKey match)
2. MetaNetX ID lookup to fetch ΔfG

Gap breakdown against MetaNetX availability:

| Gap category | Count |
|---|---:|
| Gap AND has MetaNetX ID | 3,621 (eQuilibrator cache doesn't have this MNX entry) |
| Gap AND NO MetaNetX ID | 5,018 (Andrew's script can't even attempt the lookup) |

## Approach 2 investigation — structure-based decomposition

Sam's preferred approach: for compounds not in eQuilibrator's cache, compute ΔfG from structure directly using eQuilibrator's decomposition machinery. The API for this is `equilibrator_assets.local_compound_cache.LocalCompoundCache`.

**Install performed:**
- New micromamba env `equilibrator` (Python 3.10)
- conda-forge: `equilibrator-api 0.6.0`, `equilibrator-cache 0.6.2`, `component-contribution 0.6.0`, `openbabel`, `rdkit`
- pip: `equilibrator-assets 0.6.0`
- Default compound cache: downloaded from Zenodo (1.34 GB `compounds.sqlite` + 68 MB `cc_params.npz`), stored at `~/.cache/equilibrator/` and copied to a working `compounds.sqlite`

**The blocker: cxcalc (ChemAxon Marvin) is required.**

`LocalCompoundCache.__init__` sets `self.read_only = True` unconditionally when `CHEMAXON_STATUS != 0`, and `get_compounds()` at line 171 unconditionally raises when `read_only` is true. Monkey-patching `lc.read_only = False` gets past the guard, but the underlying `create_compounds` code path still invokes `chemaxon` first (before its own "bypass" branch), and this throws `ChemAxonNotAvailableError: No such file or directory (Please ensure that Marvin cxcalc is properly installed as described at https://chemaxon.com/.)`.

The `bypass_chemaxon=True` and `specified_pkas={smi: [pkas]}` parameters exist but are **fallbacks for when cxcalc IS present** — not paths for when cxcalc is absent. They handle "cxcalc failed on this specific compound" cases; they don't replace cxcalc entirely.

Note: `chemaxon.get_atom_bag()` — despite the module name — uses OpenBabel (`pybel`) and does NOT need cxcalc. It's the pKa-decomposition step and other functions in the same module that require cxcalc.

## The ModelSEED pKa/pKb data is already there

Coverage of the per-compound pKa/pKb fields in `compound_*.json`:

- 32,223 compounds with non-empty `pka` (71% of total, ≈95% of non-R structured)
- 33,964 with non-empty `pkb`

Format:
- Top-level `pka` / `pkb`: semicolon-separated triples `"fragment_idx:atom_idx:value"`
- Per-tool dict `pkas`: keyed by predictor name (e.g., `{"Marvin": {"pKa": "...", "pKb": "..."}, "MolGpKa": {...}}`)

The `pkas` per-tool dict already has both:
- **Marvin** (ChemAxon-derived, from an earlier run when there was a license — using these values doesn't require a current license)
- **MolGpKa** (open-source alternative)

So Sam has the raw material to feed eQuilibrator's `specified_pkas` parameter — if the cxcalc block can be avoided.

## Paths forward

1. **Install ChemAxon Marvin under academic license.** ANL should qualify. Turnaround ~a day. Once installed, everything in the approach-2 plan works: `LocalCompoundCache.get_compounds([smi_list], specified_pkas={smi: [our_pkas]})` computes structure-based ΔfG for compounds not in eQuilibrator's standard cache. This is the intended path.

2. **Custom pipeline bypassing LocalCompoundCache.** Use `equilibrator_cache` directly to construct `Compound` objects from our SMILES + our MolGpKa/Marvin pKas, insert into a SQLite DB, and plug into `ComponentContribution`. Requires internalising eQuilibrator's Compound schema and decomposition API. Realistic estimate: 2–3 days of focused development.

3. **(Chosen 2026-08-04) Fall back to an extended approach 1.** Extend Andrew's `Retrieve_eQuilibrator_Compound_Energies.py` with a structure-based InChIKey fallback: for each ModelSEED non-R structured compound, first try MetaNetX ID (unchanged), then try `inchikey:<inchikey>` accession. Doesn't reach compounds truly missing from the eQuilibrator cache, but recovers the 3,621+ compounds whose InChIKeys ARE in the cache but whose MetaNetX IDs are not resolvable. Per the user's constraint: InChIKey is structure-derived (a hash of the InChI), so this satisfies the "must be structure-based" requirement.

## Decision

Path 3 for this cycle. Paths 1 and 2 stay open as follow-ups. If Sam obtains a ChemAxon Marvin academic license later, revisit Path 1 for the truly-uncached compounds.

## Environment reproducibility

To recreate the environment:

```
micromamba create -n equilibrator -c conda-forge -y \
    python=3.10 pip openbabel rdkit pandas numpy
micromamba install -n equilibrator -c conda-forge -y \
    equilibrator-api equilibrator-cache
/scratch/seaver/micromamba/envs/equilibrator/bin/pip install equilibrator-assets
```

The Zenodo compound cache (`compounds.sqlite` + `cc_params.npz`) autodownloads to `~/.cache/equilibrator/` on first `LocalCompoundCache()` instantiation. Total ~1.4 GB.
