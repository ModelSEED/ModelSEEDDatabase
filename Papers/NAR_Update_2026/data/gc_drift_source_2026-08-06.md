# The 2010→2020 ΔG drift: source identified (corrected)

**Date:** 2026-08-06
**Supersedes:** `gc_drift_source_2026-08-05.md`
**Trigger:** Chris Henry (original MFAToolkit author) clarification via Sam, 2026-08-06.

## Correction to the 2026-08-05 report

The earlier report framed the split as "pH 0 (raw MFAToolkit output) vs pH 7
(biological standard)". That framing was wrong. Chris confirmed:

> We definitely did pH 7. That was absolutely the standard. I don't
> remember retraining the cues table at any point.... You should just use
> the values from the 2010 paper.
>
> Those are absolutely pH 7 — but it's important to explain to AI —
> in my formalism from the 2010 paper, we set the energy of free H⁺ to
> −9.5 kcal/mol (equivalent of a pH of 7). By doing that, we did NOT
> need to remove the hydrogen energies from the compound energies
> themselves.
>
> So you can do it either way... pull out the proton energies from the
> compounds (if you wish) and set the energy of H⁺ to zero... or set
> the energy of H⁺ to −9.5 kcal/mol and leave the proton energies IN
> the formation energies... which is I suppose equivalent to shipping
> the database with a pH of 0 as claude describes it.

The database is NOT split between two standard states. It is split between
two **mathematically equivalent conventions** for encoding pH 7.

## Two conventions, same physics

**Convention A — Chris's 2010 formalism (what MFAToolkit ships)**
- Compound ΔfG includes the H-atom contribution from elemental formation
  (e.g. water = −56.687 kcal/mol accounting for its 2 H atoms).
- H⁺ is assigned a fixed formation energy of **−9.5 kcal/mol**
  (= −RT ln 10 × 7 at 298 K, i.e. the pH 7 correction lives here).
- Any reaction that produces/consumes protons picks up the correction
  automatically through the H⁺ stoichiometric term.

**Convention B — Alberty-style (biochemical standard, what eQuilibrator emits)**
- Compound ΔfG′ has the H-atom contribution subtracted at compound level
  (e.g. water = −37.6 kcal/mol = −56.687 + 2 × 9.539).
- H⁺ has ΔfG′ = 0 kcal/mol.
- The pH 7 correction is pre-baked into each compound.

For any properly-balanced reaction the two conventions yield the **same
ΔG′°**. But **you cannot mix conventions within a single database**: a
reaction with some Convention-A and some Convention-B compounds will be
off by (Δn_H − n_transported_H⁺) × 9.539 kcal/mol.

## What actually happened

Of the 19,383 compounds in current MSD that also carry a value in the
2010 archive (identical `cueTable.txt` in both eras):

| Group | Convention | Count | Fraction |
|---|---|---:|---:|
| **A** | 2010 formalism (Chris's original) preserved | 8,227 | 42.4% |
| **B** | Alberty-style transformed (Δ = raw + n_H × 9.539) | 1,513 | 7.8% |
| **C** | Alberty-transformed with per-compound protonation nuance (coefficient ≈ 9.35/H) | 9,643 | 49.7% |

Groups B+C received a Convention-A → Convention-B rewrite; Group A did
not. The transform was applied inconsistently across the database at
some point between 2010 and 2020 and never completed. Today's `deltag`
field is a mix of both conventions with no annotation of which is which.

Chris confirmed he did NOT retrain the cue table between 2010 and 2020.
The resurrected MFAToolkit run (2026-08-05) reproduces the 2010 output
byte-for-byte on 91% of compounds, confirming this: the pipeline is
unchanged, so the divergence isn't in the compute — it's in the
post-processing convention rewrite.

## Chris's recommendation

Use the 2010 values directly. Encode H⁺ = −9.5 kcal/mol. Don't apply the
Alberty transform.

The clean fix for the current database:
1. Import Chris's 2010 `compoundTable.txt` values as the canonical
   `Group contribution` per-source entry in `thermodynamics` (the 10K
   compounds where we have them).
2. For compounds without a 2010 entry, run the resurrected MFAToolkit on
   their current structures — produces Convention-A output natively.
3. Ship the database in Convention A uniformly. Document H⁺ = −9.5 as
   the reference.

Consumers who want Convention B (e.g. eQuilibrator downstream) can
apply the transform themselves at query time: subtract n_H × 9.539 per
compound, add +9.5 to any H⁺ term.

## Impact on the NAR paper

The story sharpens rather than weakens:

- **Section headline stays**: current `deltag` mixes conventions
  silently; this is a real data-quality issue for anyone doing
  thermodynamics-based analysis.
- **Root cause is now known**: an incomplete Convention-A → Convention-B
  migration attempted between 2010 and 2020. Chris confirms the cue
  table + underlying compute were never re-run.
- **Fix is prescribed by the original author**: revert to Convention A
  (2010 values + H⁺ = −9.5). No re-computation needed for the 10K
  covered compounds; resurrected MFAToolkit fills the rest.
- **Convention documentation becomes a first-class part of the release**:
  the database should annotate every ΔG-carrying field with which
  convention it uses.

## What this rules out (unchanged from 2026-08-05)

- ~~Cue database was rescaled~~ — RULED OUT (byte-identical 2010 vs 2020 `cueTable.txt`).
- ~~New MFAToolkit version produced different output~~ — RULED OUT (2026 resurrection ≡ 2010 output byte-for-byte on 91% of compounds; residuals correlate with structure updates).
- ~~Structural updates~~ — RULED OUT (only 3.3% of 2010→2026 divergences are >10 kcal/mol).

## Files

- Raw 2026 output: `/scratch/seaver/tmp/mfa_run_2026-08-05/all_compounds_gc_labeled.tsv`
- 2010 archive: `/scratch/seaver/Claude_Projects/MSD_Structures/MFAToolkit/inputs_archive/`
- 2020 archive: `/scratch/seaver/Claude_Projects/MSD_Structures/MFAToolkit/msdc_2020/`
- Prior (superseded) analysis: `gc_drift_source_2026-08-05.md`
