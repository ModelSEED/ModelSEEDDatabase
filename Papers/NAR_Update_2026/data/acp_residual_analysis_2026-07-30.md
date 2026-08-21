# ACP residual analysis — 2026-07-30

Analysis of the 260 open acyl-ACP name-matched candidates after the earlier bulk-add-47-safe-ACP round (commit `9819f3f`, 2026-07-02). Run against `dev` HEAD **`3c9a977`**.

## Pipeline applied

1. **Structural filter** — of the 260 open ACP name-matches, keep only compounds with `R` in the stored formula OR `*` in the SMILES: **179 survivors**, 81 dropped as false positives (mostly compounds with `null` formulas — abstract wildcard classes like `acyl-ACP(n+2)` — plus a few specific small molecules like `(S)-ACPA`).
2. **Formula-proposal** via `Scripts/Structures/Compute_ACP_Overrides.py` on the 179 survivors: **173 proposed rows** (13 skipped as `unexpected_wildcard_count_3` — compounds with 3 wildcards in their SMILES, outside the KEGG-style / MetaCyc-style pattern the script handles).
3. **Impact report** via `Scripts/Structures/Report_Formula_Change_Impact.py` on the 173 proposals:

    | Bucket | Count | Meaning |
    |---|---:|---|
    | newly_imbalanced | 53 | Currently balanced reactions that would break if all 173 applied |
    | newly_balanced | 2 | Currently broken reactions that would be fixed |
    | still_imbalanced | 33 | Already broken, unaffected |
    | still_balanced | 143 | Already OK, unaffected |

4. **Iterative safe-filter** — starting from the 173 proposals, iteratively drop any proposal implicated in a newly_imbalanced reaction, re-run, repeat until zero newly_imbalanced remain. Applied to model the co-application effect where a proposal is safe only when its reaction partners are also being updated in the same PR.

    | Iteration | Newly-imbalanced remaining | Proposals dropped | Proposals remaining |
    |---:|---:|---:|---:|
    | 1 | 53 | 55 | 118 |
    | 2 | 24 | 23 | 95 |
    | 3 | 21 | 19 | 76 |
    | 4 | 16 | 13 | 63 |
    | 5 | 14 | 12 | 51 |
    | 6 | 12 | 11 | 40 |
    | 7 | 10 | 10 | 30 |
    | 8 | 12 | 9 | 21 |
    | 9 | 9 | 7 | 14 |
    | 10 | 6 | 6 | 8 |
    | 11 | 6 | 5 | 3 |
    | 12 | 6 | 3 | 0 |
    | 13 | 0 | — | **0 (converged)** |

    Result: **0 safe proposals**. The residual ACPs cannot be bulk-corrected via the automated pipeline — every candidate cascades into breakage through co-application effects.

## Priority-scope overlap

- **0 of 173** proposals appear in any of the 9,000 v7.0 template reactions.
- **0 distinct priority reactions** would be touched by any of the proposals.

## Interpretation

The 173 residual ACP compounds are concentrated in specialty secondary metabolism pathways — polyketide synthase (PKS) chains, non-ribosomal peptide synthetase (NRPS) chains, and specific biosynthesis pathways (actinorhodin, mupirocin, mithramycin, and related). The reactions in these pathways are heavily interconnected chains where each ACP-bound intermediate is a substrate for the next reaction, so a single formula change cascades into breakage several steps downstream.

The 47-safe result from the earlier bulk-add round represented the "low-hanging fruit" — ACP compounds whose formula deltas did not propagate to their reaction partners. The 0-safe result here represents the harder residuals — compounds whose reactions form tight cycles that admit no single-compound-at-a-time fix.

## Implications for the 2026 paper

- **Priority-scope claim is unchanged:** all priority-scope acyl-ACP compounds already have pantetheine-inclusive overrides (via commit `9819f3f`). The 173 residuals are entirely outside priority scope.
- **Beyond-priority ACP curation is future work.** The Methods paragraph on ACP standardisation can honestly report: "570 acyl-ACP compounds now carry pantetheine-inclusive formula overrides. Automated extension of the framework to the 260 name-matched residuals identified 173 compounds where a formula proposal could be computed from source SMILES; iterative safe-filter analysis then found zero proposals could be applied without newly imbalancing at least one currently-balanced reaction, reflecting the tight interconnection of ACP-bound intermediates in specialty secondary-metabolism pathways (polyketide synthase chains, NRPS chains, and specific biosynthesis routes including actinorhodin, mupirocin, and mithramycin). Systematic curation of these residuals will require coupled-update reasoning that considers multiple compounds together, and is deferred to future work."
- **The 47-safe / 0-safe contrast** is itself a defensible methods observation worth a line in Results.

## Next step for §5

Pivot from ACP-residuals to the two remaining carrier classes:
- **Biotinyl-BCCP** — small class (~4-5 true carrier compounds after inspection), hand-authorable override file.
- **Lipoyl-carrier** — 82 name-matched compounds needing per-compound audit (which are free small molecules vs which are protein-bound forms).
