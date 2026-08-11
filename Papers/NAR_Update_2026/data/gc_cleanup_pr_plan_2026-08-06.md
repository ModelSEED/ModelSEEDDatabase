# GC energy cleanup PR — plan

**Date:** 2026-08-06 (revised)
**Follows:** `gc_drift_source_2026-08-06.md`
**Author:** Sam Seaver + Claude (planning session)
**PR base:** dev (per repo convention)

## Goal

Ship the ModelSEED biochemistry database with a **single, uniform,
documented convention** for `Group contribution` entries in every
compound and reaction's `thermodynamics` dict.

The target convention is Chris Henry's original (as authored in the
2010 paper and produced by MFAToolkit natively):

- Compound ΔfG *includes* full hydrogen accounting
- `ΔfG(H⁺) = −9.5 kcal/mol` (encodes the pH 7 correction)
- No per-compound Legendre transform

## Scope decision (locked 2026-08-06)

**Do NOT touch the top-level scalars `deltag`, `deltagerr`,
`reversibility`.** They retain their historical values (which are a
mix of Convention A/B via the `Promote_Reaction_Thermodynamics_to_
Canonical.py` tier policy). Retirement of the top-level scalars is a
future PR after the reorganization.

**All new work lands in the per-source `thermodynamics` dict.** The
`thermodynamics["Group contribution"]` entry becomes uniformly
Convention A after this PR; the other keys (`eQuilibrator`,
`dGPredictor`, `dGPredictor-ModelSEED`) are unchanged.

This preserves backwards compatibility (any consumer reading top-level
`deltag` today keeps working) while delivering a clean per-source
Convention A view for consumers that read `thermodynamics`.

## Approach: fresh MFAToolkit run over every compound (no import of 2010 values)

Per Sam's decision (2026-08-06): don't import the historical 2010
`compoundTable.txt` values — regenerate every compound's GC energy
from its current SMILES using the resurrected MFAToolkit. The 2010
values become a *validation reference* only, not a data source.

Rationale:
- The pipeline is deterministic and produces Convention A natively.
- Current structures have drifted from 2010 for some compounds
  (curation updates, RDKit recanonicalization). Regenerating uses the
  current structure of each compound, not the frozen 2010 snapshot.
- Sidesteps the "which compounds have 2010 values" bookkeeping — every
  structured compound gets a fresh value in the same run.
- Verified: resurrected MFAToolkit reproduces 2010 output byte-for-byte
  on 91% of compounds where structure hasn't changed
  (`mfa_run_2026-08-05/all_compounds_gc_labeled.tsv`).

## Scope

| Compound class | Count | Action |
|---|---:|---|
| Structured (SMILES present, no R) | 30,479 | Rerun MFAToolkit; write result to `thermodynamics["Group contribution"]` |
| MFAToolkit-labelable subset | 25,010 | Get numeric ΔfG |
| Labels but can't assign all atoms | 5,343 | Sentinel + note |
| Excluded: giant polymers (SMILES > 500 chars) | 126 | Sentinel; note in exclusion file |
| No structure / has R-group | 15,229 | Sentinel; unchanged |
| **Total compounds** | **45,708** | |

Reactions: `Update_Reaction_GroupContribution_Energies.py` aggregates
per-compound → per-reaction sums under `thermodynamics["Group
contribution"]`. All 56,012 reactions are re-processed; those with any
sentinel component get a sentinel result.

## Two GitHub repos, one relationship

Per Sam's decision:

- **`ModelSEED/MFAToolkit`** (new, separate repo) — tool source + cue
  database + build system. Tagged releases. Owns the C++ code, the
  Makefile, the cue metadata + `.gds` structure files, and the
  DefaultDBSpec / Defaults etc.

- **`ModelSEED/ModelSEEDDatabase`** (this repo) — commits the raw
  MFAToolkit output TSV as a data source under
  `Biochemistry/Thermodynamics/ModelSEED/`, plus a manifest pinning the
  MFAToolkit version (git SHA or release tag) used to produce it. The
  `Update_Compound_GroupContribution_Energies.py` script reads that
  TSV to populate compound `thermodynamics`.

Reproducibility guarantee: anyone with the DB can (a) see which
MFAToolkit release produced the current values, (b) rerun that release
against the current compound SMILES, and (c) diff.

## Phased plan (each phase is one commit; PR bundles them)

### Phase 1 — Convention documentation

**Files added:**
- `Biochemistry/Thermodynamics/README.md` — document Convention A
  explicitly. Include the H⁺ = −9.5 kcal/mol rule, the ΔfG-includes-H
  rule, worked reaction-ΔG example showing Convention A vs Convention B
  agree, and the "top-level scalars are historical, `thermodynamics`
  dict is current" boundary.

**Why first:** the rest of the PR only makes sense in the context of
this doc. Also lets reviewers push back on framing before regeneration.

### Phase 2 — MFAToolkit separation + output TSV storage

**In `ModelSEEDDatabase` (this PR):**
- `Biochemistry/Thermodynamics/ModelSEED/MFAToolkit_version.txt` —
  one-line manifest with the MFAToolkit git SHA / release tag pinned
- `Biochemistry/Thermodynamics/ModelSEED/<something>_MolAnalysis.tbl` —
  raw MFAToolkit compound output (the TSV consumed by
  `Update_Compound_GroupContribution_Energies.py`); ~1 MB
- No changes to the tool location within this repo (tool is already
  separate at `/scratch/…/MFAToolkit/` and will move to the new
  `ModelSEED/MFAToolkit` repo as part of this PR's release notes)

**In `ModelSEED/MFAToolkit` (separate repo, out of PR scope):**
- Tool source, Makefile, cue db, DefaultDBSpec — everything currently
  under `MFAToolkit/`
- First tagged release: `v2.0.0-post-resurrection` (with H2+H3 fix,
  DefaultDBSpec, current cue database bundled)

### Phase 3 — Regenerate compound GC energies via fresh MFAToolkit run

Run the tagged MFAToolkit release against all 30,479 structured
compounds. Write results into the MolAnalysis TSV committed in Phase 2.

Then run `Update_Compound_GroupContribution_Energies.py` — reads the
TSV, populates `thermodynamics["Group contribution"]` on every
compound. Also set `cpd00067` (H⁺) `thermodynamics["Group contribution"]
= [-9.5, 0.0]` explicitly (the load-bearing constant of Convention A;
`deltagerr = 0` because it's a defined reference, not an estimate).

**Does NOT touch:** top-level `deltag` / `deltagerr` on compounds.
Those keep their historical (mixed-convention) values until the future
retirement PR.

**Estimated diff:** every structured compound's `thermodynamics[
"Group contribution"]` entry rewrites. Many by <5 kcal/mol from
current dev value (Convention A stayed put), some by ~+n_H × 9.539
(where dev's GC entry had accidentally been Legendre-transformed).

### Phase 4 — Regenerate reaction GC energies + per-source operators

Rerun `Update_Reaction_GroupContribution_Energies.py` — sums per-
compound GC energies into per-reaction. Because compound GC values are
now uniformly Convention A, reaction GC values are too.

Rerun `Add_Reaction_Thermodynamics_Operators.py` — refreshes the
per-source `[dg, dge, operator]` triples on all sources (the operators
were the third element added by PR #263; this pass recomputes them
from the current `[dg, dge]` using the shared heuristic).

**Does NOT run:** `Estimate_Reaction_Reversibility.py` (writes
top-level `reversibility`; we're not touching those). Does NOT run:
`Promote_Reaction_Thermodynamics_to_Canonical.py` (writes top-level
`deltag`/`deltagerr`; also not touching).

**Effect on the per-source operators:** Since compound-level GC values
changed for many compounds, the per-reaction GC dg values change, which
means the per-reaction GC operators are recomputed and may flip for
some reactions (~1-2K expected, matching the H2+H3 fix's known impact
plus the Convention A ⇄ B compound-side corrections).

### Phase 5 — Test suite + verification

**`Scripts/Tests/test_gc_convention.py`** — new script. Picks ~10-15
landmark compounds (water, ATP, ADP, glucose, pyruvate, NAD, NADH,
Pi, PPi, NH3, CO2, H+, ...) and asserts:
- `thermodynamics["Group contribution"]` matches the historical 2010
  archive value within ±0.05 kcal/mol
- `cpd00067` `thermodynamics["Group contribution"][0]` == −9.5
- Formula-derived hydrogen count matches — a sanity check that
  Convention A structure semantics are preserved

Add to CI. Guards against future accidental convention drift.

**`Scripts/Tests/test_reaction_direction.py`** — no change to the code;
after Phase 4 lands, its baseline needs bumping. When the PR merges,
consumers running this script against the new `origin/dev` will pass.

### Phase 6 — Paper writeup update

Update `Papers/NAR_Update_2026/PAPER_2026_PLAN.md` and the drift-source
report to reflect the completed cleanup. The "before" state (mixed
conventions in `thermodynamics["Group contribution"]`) becomes the
paper's data-quality motivation; the "after" state (uniform Convention
A + documented + reproducible via bundled MFAToolkit) becomes the
contribution.

Note that top-level scalars remain historical/mixed; call this out
explicitly as scheduled for a future release.

## PR structure

**Single PR with phased commits** as above. Each commit self-contained
and independently reviewable. PR description leads with the drift-source
report + Chris's guidance, so reviewers see the motivation.

**Target branch:** dev
**PR title:** `Uniform Convention A in thermodynamics["Group contribution"]; MFAToolkit re-run + convention doc`

## Verification checklist (before opening the PR)

1. Landmark compounds match historical 2010 within ±0.05 kcal/mol:
   water, ATP, ADP, glucose, pyruvate, NAD, NADH — in
   `thermodynamics["Group contribution"]` (not top-level `deltag`)
2. `cpd00067` (H⁺) has `thermodynamics["Group contribution"][0] = -9.5`
3. Rerunning the pipeline is idempotent (byte-identical on second run)
4. Top-level `deltag` / `deltagerr` / `reversibility` are BYTE-IDENTICAL
   to pre-PR values on every compound and reaction (should not have
   changed at all)
5. Sample 20 random reactions: `thermodynamics["Group contribution"]` is
   internally consistent (sum of per-compound GC values with signs
   matches the reaction entry within rounding)
6. Numeric-GC coverage count before/after ≈ similar (~25K compounds
   with numeric values)

## Files that end up committed

New:
- `Biochemistry/Thermodynamics/README.md` — convention doc
- `Biochemistry/Thermodynamics/ModelSEED/MFAToolkit_version.txt` —
  pinned version manifest
- `Biochemistry/Thermodynamics/ModelSEED/<name>_MolAnalysis.tbl` — raw
  MFAToolkit output
- `Scripts/Tests/test_gc_convention.py` — landmark test

Modified:
- Every `Biochemistry/compound_*.json` — 30K entries with new
  `thermodynamics["Group contribution"]` values
- Every `Biochemistry/reaction_*.json` — regenerated
  `thermodynamics["Group contribution"]` values + refreshed operators
  in the full per-source thermodynamics dict
- `Biochemistry/reaction_*.tsv` — regenerated
- `Scripts/Thermodynamics/README.md` — cross-link to convention doc,
  update pipeline description

Unmodified (deliberately, backwards-compat):
- Top-level `deltag`, `deltagerr`, `reversibility` on every compound
  and reaction
- `Promote_Reaction_Thermodynamics_to_Canonical.py` (still lives in the
  repo for anyone wanting to regenerate the historical scalars)
- Any other per-source thermodynamics keys (`eQuilibrator`, `dGPredictor`,
  `dGPredictor-ModelSEED`) — those keep their current values

## Remaining open questions

1. **Draft PR vs open PR** — draft first for reviewer preview, then
   promote when Phase 5 tests pass? Recommend draft first.
2. **Deprecation notice for top-level scalars** — do we add a comment
   or doc line saying "top-level `deltag` is scheduled for retirement
   in a future release, prefer `thermodynamics[<source>][0]`"?
3. **MFAToolkit repo name** — `ModelSEED/MFAToolkit` (matches the
   existing tool name; simplest) vs `ModelSEED/mfatoolkit-cxx` (marks
   the C++ implementation explicitly, in case a Python rewrite is
   ever wanted). Recommend `ModelSEED/MFAToolkit`.

## Follow-up PRs (out of scope for this one)

1. **Top-level scalar retirement.** Once consumers have migrated to
   reading `thermodynamics[<source>]`, remove the top-level
   `deltag`/`deltagerr`/`reversibility` fields (or redefine them as
   dicts).
2. **Convert non-GC entries to Convention A.** eQuilibrator and
   dGPredictor sources emit Convention B. A follow-up could apply
   `n_H × 9.539` per-compound to bring them to Convention A too,
   yielding a fully-uniform per-source dict. Not this PR because it
   touches many more values with per-compound bookkeeping.
3. **Deprecate `Promote_Reaction_Thermodynamics_to_Canonical.py`.**
   Once top-level scalars retire.
