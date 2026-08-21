# GC ΔG divergence: 2010 archive vs. current MSD

**Date:** 2026-08-05
**2010 source:** `MFAToolkit/inputs_archive/.../compoundTable.txt` (16,274 rows, 10,164 with real ΔG)
**Current source:** `ModelSEEDDatabase/Biochemistry/compound_*.json` (`deltag` field; all entries annotated `notes=["GC","EQ","EQU"]`)

## Headline

Of the **9,306** compounds where both vintages have a numeric ΔG (formulas match for 7,123, differ for 2,183):

- **Median |shift| = 48.7 kcal/mol; mean = 101.7 kcal/mol; 95th percentile = 393 kcal/mol.**
- Distribution is strongly **bimodal**: **3,886 compounds (42%) are essentially unchanged** (|Δ| < 0.5 kcal/mol; 3,701 are byte-identical), and **5,373 (58%) shifted by > 5 kcal/mol** — very little in between.
- The direction is systematic: **1,981 sign flips are neg → pos, only 9 are pos → neg**.
- **Formula was UNCHANGED** for **73% (3,899 / 5,373)** of the large-divergence set — GC-method drift dominates, not structural updates.

## Coverage (all 16,274 shared IDs)

| Bucket                                     |     N |
| ------------------------------------------ | ----: |
| Current numeric ΔG present                 | 10,321|
| Current sentinel / missing                 |  5,953|
| Formula match (both non-empty)             | 10,934|
| Formula changed                            |  4,308|
| Current formula missing                    |    935|
| Both dg present (analysis set)             |  9,306|

## Distribution of |current − 2010|

```
[0,    0.5) kcal  3886  ############################################################
[0.5,  1)      8
[1,    2)     16
[2,    5)     23
[5,   10)    272  ####
[10,  20)     88  #
[20,  50)    395  ######
[50, 100)   1342  ####################
[100,500)   2988  ##############################################
[500,10k)    288  ####
```

Quantiles of |diff|: q25=0.00, q50=48.7, q75=137.1, q90=267.5, q95=393.2, q99=686.6, max=3007 kcal/mol.
Signed: **57.7% current > 2010** (mean +101.1, median +48.5) vs. 2.5% current < 2010; the remainder are the 42% unchanged tail. The 2020 GC recomputation is monotonically less negative / more positive than the 2010 one for anything it actually touched.

## Small-molecule sanity checks

| cpd      | name       | 2010    | current | diff    |
| -------- | ---------- | ------: | ------: | ------: |
| cpd00007 | O2         |   +3.92 |   +3.92 |   0.00  |
| cpd00011 | CO2        |  −92.26 |  −92.26 |   0.00  |
| cpd00067 | H+         |   −9.53 |    0.00 |  +9.53  |
| cpd00001 | H2O        |  −56.69 |  −37.54 | +19.15  |
| cpd00009 | Phosphate  | −261.97 | −252.51 |  +9.46  |
| cpd00013 | NH3        |  −18.97 |  +19.05 | +38.02  |
| cpd00033 | Glycine    |  −87.73 |  −41.92 | +45.81  |
| cpd00027 | D-Glucose  | −218.28 | −102.12 |+116.16  |
| cpd00002 | ATP        | −673.85 | −548.85 |+125.00  |
| cpd00008 | ADP        | −465.85 | −340.04 |+125.81  |

Even for the most basic metabolites the current values sit ~10–125 kcal/mol above 2010. The uniform positive-signed drift and the tight ATP ≈ ADP shift (+125.0 vs +125.8 kcal/mol) is highly suggestive of a re-referenced cue database (element-formation-energy convention change) rather than random per-compound recuration.

## Sign flips (n = 1,990; 1,673 formula-unchanged)

Direction is essentially one-way (**1,981 neg→pos, 9 pos→neg**). Small-molecule cases where the qualitative direction of stability inverted:

| cpd      | name                        | 2010    | current | formula same? |
| -------- | --------------------------- | ------: | ------: | :-----------: |
| cpd00013 | NH3                         |  −18.97 |  +19.05 |      yes      |
| cpd00239 | H2S                         |   −6.66 |  +12.31 |      no       |
| cpd01530 | Dichloromethane             |  −15.80 |   +3.68 |      yes      |
| cpd03957 | Ethylene oxide              |  −15.84 |   +9.41 |      yes      |
| cpd04749 | Pteroic acid                |  +11.37 |   −0.82 |      yes      |
| cpd02506 | 5-Amino-4-imidazoleCA       |   +0.09 |   −5.28 |      yes      |

## Top-10 large |diff| with formula UNCHANGED (real GC drift, not structural update)

| cpd      | name (trunc.)                            | 2010     | current  | |diff|   | formula              |
| -------- | ---------------------------------------- | -------: | -------: | -------: | -------------------- |
| cpd15519 | (O16 antigen)x4 undecaprenyl-diphosphate | −3307.50 |  −300.68 |  3006.82 | C191H310N4O107P2     |
| cpd15517 | (O16 antigen)x3 undecaprenyl-diphosphate | −2519.77 |   −49.43 |  2470.34 | C157H255N3O82P2      |
| cpd15516 | (O16 antigen)x2 undecaprenyl-diphosphate | −1732.04 |  +201.81 |  1933.85 | C123H200N2O57P2      |
| cpd15456 | (enterobacterial common antigen)x3 uPP   | −1553.02 |  +363.45 |  1916.47 | C127H198N9O52P2      |
| cpd15491 | KDO(2)-lipid IV(A) w/ palmitoleoyl       | −1403.83 |  +349.25 |  1753.08 | C100H178N2O38P2      |
| cpd15455 | (enterobacterial common antigen)x2 uPP   | −1087.54 |  +477.45 |  1564.99 | C103H162N6O37P2      |
| cpd03492 | Undecaprenyl-diphospho-N-acetylmuramoyl… |  −709.35 |  +778.20 |  1487.55 | C94H155N9O25P2       |
| cpd03491 | Undecaprenyl-diphospho-N-acetylmuramoyl… |  −751.79 |  +718.00 |  1469.79 | C94H153N8O26P2       |
| cpd03495 | Undecaprenyl-diphospho-N-acetylmuramoyl… |  −831.43 |  +629.27 |  1460.70 | C95H152N8O28P2       |
| cpd09661 | alpha-Semegma mycolic acid               |   +39.62 | +1491.03 |  1451.41 | C77H149O3            |

## Formula-changed vs formula-unchanged (both dg present, |diff|>5)

|                    | total | formula changed | formula unchanged |
| ------------------ | ----: | --------------: | ----------------: |
| |diff| > 5         | 5,373 |           1,415 |             3,958 |
| |diff| > 10        | 5,101 |           1,202 |             3,899 |

Median |diff| is similar for both subsets (51 kcal/mol same-formula vs 37 kcal/mol formula-changed). Conclusion: **structural updates are NOT the dominant explanation** — the recalibration touches same-formula compounds just as hard.

## Structural family of the same-formula divergent set (|d|>10, N = 3,899)

Rough classification by formula string (multiple tags allowed):

| pattern (from formula)             |    N |
| ---------------------------------- | ---: |
| long C-chain (C ≥ 10)              | 2,459|
| contains N                         | 1,614|
| contains P                         |   190|
| contains S                         |   119|
| contains a metal (Fe/Cu/Mn/…)      |     0|
| plain (no P/S/N, C<10)             |   752|

Large-divergence, same-formula compounds are heavily enriched in **long carbon chains** — polyprenyl carriers (undecaprenyl-*), lipopolysaccharides (KDO2-lipid A, O-antigens), mycolic acids, and carotenoids (Thermobiszeaxanthin, Physalien, Helenien). These are exactly the structural classes where the group-contribution decomposition sums over many repeat units, so a small per-cue rescaling accumulates into ~1000+ kcal/mol shifts.

## Bottom line

The 2010 → 2020 shift is **not small and boring** — it is large, systematic, and one-directional. About 42% of shared compounds carried the exact 2010 value forward untouched; the other ~58% were recomputed and moved by tens to hundreds of kcal/mol, almost always toward less-negative ΔG. Even ATP, ADP, and glucose moved +100–130 kcal/mol without any formula change. The pattern (uniform-magnitude shifts on chemically related pairs like ATP/ADP, one-way sign flips, amplification on long-chain compounds) is the fingerprint of a **rescaled / re-referenced cue table**, not per-compound recuration and not structural updates.
