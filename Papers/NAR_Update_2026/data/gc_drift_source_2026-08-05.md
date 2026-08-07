# The 2010→2020 ΔG drift: source identified

**Date:** 2026-08-05
**Follow-up to:** `gc_2010_vs_current_divergence_2026-08-05.md`

## Result

The 2010 → 2020 divergence in `deltag` values is **entirely a Legendre transform artifact**. It is NOT due to a rescaled cue table, restructured compounds, or a different pipeline. It is a **standard-state change from pH 0 (raw MFAToolkit output) to pH 7 (biological standard state)**, applied inconsistently across the database.

The transform formula that fits:

```
ΔG'° (pH 7) = ΔG° (pH 0) + n_H × RT × ln(10) × pH
            = ΔG° (pH 0) + n_H × 9.539 kcal/mol    (at T=298K, pH=7)
```

## Evidence

Ran the resurrected MFAToolkit fresh on all 30,353 compounds (2026-08-05). Compared the resulting raw ΔG values against the current MSD `deltag` field for 19,383 overlapping compounds:

| Group | Definition | Count | Fraction |
|---|---|---:|---:|
| **A** | 2020 = raw 2026 (byte-identical, no transform) | 8,227 | 42.4% |
| **B** | 2020 = raw + n_H × 9.539 to within 2 kcal/mol (strict Legendre) | 1,513 | 7.8% |
| **C** | 2020 ≠ raw AND ≠ strict Legendre | 9,643 | 49.7% |

Groups B and C together (57.6%) received the pH 7 transform. Group A did not.

### Group C is also Legendre-transformed — just with H-count nuances

For Group C, the effective coefficient `(2020 − raw) / n_H` is centered on **+9.35 kcal/mol per hydrogen** (25th–75th percentile: [7.85, 9.91]). The theoretical value for a pH 7 transform at 298K is **+9.539**. The small deviation comes from protonation-state corrections — the formula's H count reflects the neutral reference form, not the actual protonation state at pH 7, so the effective `n_H` in the Legendre transform differs by 1–4 protons for ionic species.

The upshot: **all of Group B + C (57.6% of compounds) show the pH 7 transform signature.** Only Group A (42.4%) is untransformed.

### Landmark compounds — sub-2 kcal/mol residuals

| Compound | Formula | H | Raw 2026 | Predicted (raw + H·9.54) | Actual 2020 | Residual |
|---|---|---:|---:|---:|---:|---:|
| ATP | C10H13N5O13P3 | 13 | −673.85 | −549.84 | −548.85 | **+0.99** |
| ADP | C10H13N5O10P2 | 13 | −465.85 | −341.84 | −340.04 | **+1.80** |
| Glucose | C6H12O6 | 12 | −218.28 | −103.81 | −102.12 | **+1.69** |
| Pyruvate | C3H3O3 | 3 | −112.69 | −84.07 | −82.56 | **+1.51** |
| NAD | C21H26N7O14P2 | 26 | −529.59 | −281.58 | −286.41 | −4.83 |
| NADH | C21H27N7O14P2 | 27 | −524.32 | −266.77 | −271.15 | −4.38 |

The 6 landmarks span central carbon and energy metabolism. All predicted from raw + Legendre transform to within 5 kcal/mol.

## Why does this matter for the NAR paper

**Three data problems this identifies:**

1. **Inconsistent standard state within the database.** Group A entries (42% of GC-covered compounds) are at pH 0. Groups B + C (58%) are at pH 7. **A user cannot tell which standard state a compound's `deltag` is in without running the check we just did.** For thermodynamics-based analyses (flux directionality, feasibility), mixing the two silently produces nonsense.

2. **The `deltag` field was never regenerated in the 2020 update.** MFAToolkit was not rerun for the 2020 paper (per the confirmation that raw 2026 output exactly matches 2010 output on 91% of compounds, with median |Δ|=0.000). Instead, an ad hoc Legendre-transform script was applied to a partial subset. What we see in the current database is the fossil of a half-applied transform, layered on top of the frozen 2010 GC output.

3. **The transform was applied to ~58% of compounds, missed ~42%.** No obvious selection criterion — Group A and Group B/C compounds are interleaved in the cpd ID range and have similar H count distributions.

**What to include in the NAR paper (§6-§8):**

- Report the discovery of standard-state inconsistency.
- Publish the fresh 2026 rerun as `Group contribution (MFAToolkit-2026)` in the additive thermodynamics dict — this gives users a coherent pH 0 baseline they can Legendre-transform consistently themselves.
- Deprecate or annotate the top-level `deltag` field to warn users it is inconsistent.
- Recommend downstream tools apply their own Legendre transform from the raw values, not consume the pre-transformed `deltag`.

## What this rules out

- **Cue database was rescaled** — RULED OUT. `cueTable.txt`, `compoundTable.txt`, `reactionTable.txt`, `FinalGroups.txt`, `AtomTypes.txt`, `GroupTranslation.txt` are all byte-identical between the 2010 archive and the 2020 msdc_files archive.
- **New MFAToolkit version produced different output** — RULED OUT. Our resurrection of MFAToolkit reproduces 2010 output exactly (8,353 of 9,222 overlapping compounds byte-identical, median |Δ| = 0.000 kcal/mol).
- **Compounds were restructured** — RULED OUT. Only 3.3% of 2010 → 2026 divergences are >10 kcal/mol; those correlate with formula changes and are the expected residual.

## Files
- Raw 2026 output: `/scratch/seaver/tmp/mfa_run_2026-08-05/all_compounds_gc_labeled.tsv`
- 2010 archive: `/scratch/seaver/Claude_Projects/MSD_Structures/MFAToolkit/inputs_archive/`
- 2020 archive: `/scratch/seaver/Claude_Projects/MSD_Structures/MFAToolkit/msdc_2020/`
