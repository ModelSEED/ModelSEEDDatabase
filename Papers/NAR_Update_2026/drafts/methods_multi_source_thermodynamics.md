# Methods — Multi-source thermodynamics (draft)

**Target section:** Materials and Methods, "Multi-source thermodynamics" (see `PAPER_2026_SKELETON.md`).
**Guide reference:** `PAPER_2026_GUIDE.md` §6.
**Status:** GC subsection filled in against the committed rebuild; eQuilibrator / dGpredictor / OpenTECR versions still `[TBD]`; direction-heuristic rule set still `[TBD]`.

---

The 2020 release combined two thermodynamics sources: eQuilibrator as the primary estimator of Δ<sub>f</sub>G′ and Δ<sub>r</sub>G′, and a Group Contribution (GC) method as a fallback for compounds and reactions eQuilibrator could not estimate. For this update we integrate four sources, each providing an independent estimate wherever its coverage permits, and store per-source `[energy, error, operator]` triples on every compound and reaction — no source is silently overwritten. The top-level `deltag` / `deltagerr` / `reversibility` scalars are still populated from a tier + lowest-uncertainty policy, but are being deprecated in favor of the per-source records.

**eQuilibrator refresh.** [TBD version] was run against the 2026 compound set at the reported conditions below. Coverage: [TBD].

**Reported conditions, and why they are not either eQuilibrator preset.** All eQuilibrator-derived values in this release are reported at **pH 7.0, pMg 3.0 (1 mM free Mg²⁺), ionic strength 0.25 M and 298.15 K**. This is a deliberate combination and matches neither of the presets the software ships:

| | pH | pMg | I | T |
|---|---|---|---|---|
| eQuilibrator *standard* | 7.0 | 10 | 0.25 M | 298.15 K |
| eQuilibrator *physiological* | 7.5 | 3.0 | 0.25 M | 298.15 K |
| **This release** | **7.0** | **3.0** | 0.25 M | 298.15 K |

We adopt **pMg 3.0** on biological grounds. Free intracellular Mg²⁺ is of order 0.5–2 mM across cell types, and the standard preset's pMg 10 corresponds to essentially no free magnesium — a state no cell occupies, and a poor reference for energies whose purpose is to constrain metabolic models, given that ATP is predominantly Mg-bound *in vivo*. We retain **pH 7.0** on grounds of convention and continuity: it is the standard biochemical reference state, it is the basis of the Convention A formalism used for the Group Contribution values (above), and it is what the 2020 release reported. Moving to the physiological preset's pH 7.5 would shift every value and break comparability both with our own prior release and with the GC source, for no gain beyond internal tidiness.

We note for transparency that this combination was originally *inherited* rather than chosen: pre-2023 retrieval scripts set `p_h` explicitly and never set `p_mg`, so the constructor's physiological magnesium default persisted alongside a standard-preset pH. The audit described next was prompted by that discovery; the conditions are unchanged, but they are now a decision rather than an accident.

*A documented inconsistency between training and prediction.* Component contribution is fitted to the NIST TECRDB corpus, in which **2,944 of 4,456 measurements (66%) record no magnesium concentration** and the loader substitutes pMg 14 — again, effectively magnesium-free. The 1,512 measurements that *do* record it have a median pMg of 2.4 (≈4 mM). Since enzyme assays for ATP-dependent chemistry require magnesium, an unrecorded value is far more likely to mean "present and unremarkable" than "absent". The model is therefore fitted largely as though the assays were magnesium-free and queried at 1 mM.

We tested whether this matters. Three training variants were cross-validated over ten folds with an identical held-out test set in every arm, differing only in the value assumed for the 2,944 unrecorded rows: pMg 14 (status quo), pMg 2.4 (the median of rows that measured it), and pMg 3.0 (making training and prediction mutually consistent).

| training assumption | RMSE | MAE | median \|e\| | vs status quo |
|---|---|---|---|---|
| pMg 14 | 3.7636 | 2.0294 | 1.1217 | — |
| pMg 2.4 | 3.8190 | 2.0480 | 1.1412 | *p* = 0.28 |
| pMg 3.0 | 3.7666 | 2.0196 | 1.1448 | *p* = 0.89 |

Neither alternative is distinguishable from the status quo. We therefore retain the loader default and document the inconsistency rather than re-fitting: a re-fit would perturb every value in the database for no measurable improvement in held-out prediction.

The residual effect is bounded and small. Recomputing Δ<sub>r</sub>G′ at pMg 14, 3.0 and 2.4 for reactions that involve at least one Mg²⁺-binding compound (511 such compounds in the cache; **34,359 of 56,012 reactions, 61%**, touch one) gives a median difference between the training assumption and the reporting condition of **0.044 kcal/mol**, a 90th percentile of 0.98 and a maximum of 2.51. Reactions with no Mg²⁺-binding participant are unaffected to machine precision, which serves as the control. The effect is thus within the reported uncertainty for the large majority of reactions, but is not negligible for reactions whose direction call sits near the reversibility threshold.

*Correction to a prior internal result.* An earlier experiment (`--mg-from-metadata`) was recorded as showing that deriving pMg from measurement metadata degraded prediction, and was cited as settling this question. It did not: that variant derives pMg from openTECR's `Cofactor` field, and the base TECRDB has no such column, so it could only ever have modified the ~178 openTECR-sourced rows — and it was evaluated inside an augmented training table that was independently found to be harmful. It is silent on the default applied to the 2,944 original rows, which is the experiment reported above.

**Group Contribution rebuilt from MFAToolkit under Chris Henry's 2010 Convention A.** The 2020 fallback GC method had drifted from its original conventions in a subtle but load-bearing way, and the update rebuilds it from the source implementation with the drift closed. Both the drift finding and the fix are worth documenting here because they change reaction Δ<sub>r</sub>G′ estimates by up to 9.5 kcal/mol on any reaction whose two source values crossed the drift line.

*The drift finding.* Audit against the 2010 GC values (Jankowski et al., <a href="https://doi.org/10.1529/biophysj.107.124784">10.1529/biophysj.107.124784</a>) revealed that the top-level `deltag` scalar in the 2020 database was a 42/58 mixture of two mathematically equivalent but incompatible conventions for encoding pH 7. In **Convention A** (Chris Henry's original 2010 formalism, and what MFAToolkit natively produces), the free-proton formation energy is fixed at Δ<sub>f</sub>G(H⁺) = −9.5 kcal/mol (numerically = −RT ln(10) × pH at pH 7) and compound Δ<sub>f</sub>G values include their full hydrogen accounting — water is −56.687 kcal/mol. In **Convention B** (Alberty's transformed formalism, and what eQuilibrator emits), the free-proton contribution is absorbed into each compound (Δ<sub>f</sub>G(H⁺) ≡ 0) — water is instead ~−37.6 kcal/mol. Both conventions produce identical Δ<sub>r</sub>G′ for a properly-balanced reaction, but mixing them within a single reaction breaks the accounting by (Δn<sub>H</sub> − n<sub>transported H⁺</sub>) × 9.539 kcal/mol.

*The rebuild.* We resurrected the MFAToolkit source tree (build recipe and cue-database provenance pinned in `Biochemistry/Thermodynamics/ModelSEED/MFAToolkit_version.txt`), regenerated the GC values from Chris's original per-source mol-file corpora (KEGG and MetaCyc, each in Original and Charged variants), and stamped them into the database uniformly under Convention A. A per-compound resolver now takes the mean Δ<sub>f</sub>G across curated-alias entries rather than the previous lowest-pick, with an error term inflated to `sqrt(mean(σᵢ²) + var(dgᵢ))` so alias disagreement widens the reported uncertainty rather than being silently collapsed. Nine load-bearing small-molecule compounds whose Δ<sub>f</sub>G values live only in the MFAToolkit cue database — H⁺, H₂O, NH₄⁺, CO₂, HCO₃⁻, H₂O₂, H₂, O₂, H₂S — are injected at the compound level under a documented anchor table.

*Cross-validation.* Every anchor value was independently checked against eQuilibrator's Δ<sub>f</sub>G′ shifted to Convention A via the Legendre transform Δ<sub>f</sub>G<sub>A</sub> = Δ<sub>f</sub>G<sub>B</sub> − n<sub>H</sub> × 9.539. Eight of the nine agree between the two independent methods to within 0.5 kcal/mol; the ninth (H₂S) shows a 2.88 kcal/mol residual attributable to eQuilibrator's pseudo-species treatment (a pKa-weighted mix of H₂S and HS⁻ at pH 7) vs the neutral-form cue value.

*Coverage.* [TBD compounds with numeric GC Δ<sub>f</sub>G; TBD reactions with numeric GC Δ<sub>r</sub>G′] — to be filled from the current DB snapshot.

**dGpredictor retrained on ModelSEED.** dGpredictor (an ML-based Δ<sub>f</sub>G estimator) has been retrained using the ModelSEED-curated compound structures as training data. This produces a ModelSEED-native estimator rather than one trained on a competing biochemistry database, and lets us report Δ<sub>f</sub>G estimates for compounds the other three methods cannot cover. Training procedure and hyperparameters: [TBD]. Coverage: [TBD]. A quinone/quinol audit is outstanding: dGpredictor reportedly carries extreme error on these pairs (see `data/quinone_quinol_investigation_[TBD].md`), and the paper will note whether audit results warrant a correction pass or a per-class error flag.

**OpenTECR integration.** OpenTECR is a community effort to standardize experimentally-measured thermodynamic constants; we integrate its values where a mapping between OpenTECR reaction identifiers and ModelSEED reaction IDs exists. Coverage: [TBD].

**Reported uncertainties are not on a common scale.** Each source reports a σ, but the three quantities answer different questions and are calibrated differently. Group Contribution reports the propagated uncertainty on the Jankowski 2008 group energies, summed in quadrature over a molecule's groups. eQuilibrator reports the fitted covariance of the component-contribution model, ‖σ<sub>fin</sub>‖ + RMSE<sub>inf</sub>·‖σ<sub>inf</sub>‖, in which the second term is an explicit sentinel rather than an error bar. dGpredictor reports a BayesianRidge posterior predictive standard deviation. A σ of 10 kcal/mol therefore does not mean the same thing in the three columns.

We measured each against held-out observations. Group Contribution never trains on TECRDB, so measured equilibrium constants are out-of-sample for it by construction; eQuilibrator and dGpredictor were assessed by ten-fold cross-validation, refitting each fold and recording the σ the fitted model itself reports for the held-out reactions. Comparisons are restricted to reactions with no net bound-hydrogen change, where the Convention A / B difference vanishes.

| source | n | median σ | median \|error\| | median \|z\| | frac \|z\|<1 |
|---|---|---|---|---|---|
| Group contribution | 1,685 | 2.91 | 1.03 | 0.32 | 73.6% |
| dGpredictor | 3,861 | 1.51 | 0.45 | 0.29 | 87.9% |
| eQuilibrator | 3,663 | 0.10 | 0.27 | **2.27** | **28.7%** |
| *calibrated reference* | | | | *0.67* | *68%* |

Group Contribution and dGpredictor are **over-conservative** by roughly a factor of two: their error bars are wider than their errors. eQuilibrator is **over-confident** by roughly a factor of three — it reports a median σ of 0.10 kcal/mol on held-out reactions whose median error is 0.27, and only 29% of predictions fall within one stated σ where 68% should.

Two points follow. First, eQuilibrator remains the most accurate source by a clear margin (median held-out error 0.27 kcal/mol against 0.45 and 1.03); the finding concerns its stated confidence, not its estimates. Second, and more consequentially, the policy that populates the canonical `deltag` selects the **lowest-uncertainty** estimate within a quality tier — so it has been systematically preferring the source that most understates its uncertainty. Users comparing σ across sources, or relying on that selection, should treat the three scales as distinct.

We therefore do **not** rescale the stored values. The miscalibration is not a constant: for Group Contribution it varies from ×1.17 in the smallest σ quartile to ×0.27 in the largest, and is not even monotonic in σ, so no affine transform puts the three on a shared axis. Instead we publish a per-source empirical scale — the observed median and 90th-percentile error for each band of reported σ — alongside the unmodified values (`data/sigma_calibration_{gc,dg,eq}.tsv`). This preserves comparability with published eQuilibrator numbers while making the meaning of each σ recoverable.

**Direction assignment.** For each reaction with an accepted Δ<sub>r</sub>G′ (from any source), we apply the reversibility cascade documented in `Scripts/Thermodynamics/reversibility_heuristics.py` (based on Jankowski et al. and refined in the current update to handle the widened stoichiometry patterns from the H₂/H₃-carrier standardization). Output is one of `=` (reversible), `>` (irreversible left-to-right), `<` (irreversible right-to-left), or `?` (indeterminate — |Δ<sub>r</sub>G′| < error). On top of this we add **per-reaction-class heuristic overlays** — hard rules for classes where the standard thermodynamic assignment is known to disagree with textbook biology (hydrolases; cofactor-loading reactions; proton-gradient transport; electron-carrier-coupled reactions where the carrier's Δ<sub>f</sub>G is indeterminate). Rule set: [TBD; discussed but not yet codified as a single overlay file].

**Combined direction ledger.** For every reaction, the pipeline records: the Δ<sub>r</sub>G′ estimate from each of the four sources (with per-source uncertainty), the base direction assignment from each source, the heuristic overlay if applicable, and the final direction. This ledger is the input to the Results section's empirical study of how the direction-source choice affects metabolic-model output.

**Reproducibility.** MolAnalysis tables from the MFAToolkit runs (four total: KEGG × MetaCyc × Original × Charged) are committed under `Biochemistry/Thermodynamics/ModelSEED/`. A landmark regression test (`Scripts/Tests/test_gc_convention.py`) guards the 9 anchor values, 12 central-metabolism compound landmarks (ATP, NAD, glucose, …), and 2 reaction spot-checks against the water canary that would catch a silent flip back to Convention B.

---

## Open loose ends flagged during drafting

- **The three sources report at three different condition frames**, and the paper should say so plainly before any head-to-head comparison. GC is Convention A at pH 7 with no magnesium term at all; eQuilibrator is Alberty-transformed at pH 7.0 / pMg 3.0; dGpredictor applies no Legendre transform whatsoever, so its values sit at whatever mixture of pH and ionic strength dominates its training rows. Compound Δ<sub>f</sub>G are consequently *not* comparable across GC and the other two without the Convention A↔B transform, and even reaction Δ<sub>r</sub>G′ agree only where the explicit proton stoichiometry matches the bound-hydrogen change. This is a larger comparability problem than the pMg question resolved above.

- Version numbers for eQuilibrator, dGpredictor, and OpenTECR need to be captured for reproducibility. Add to the Methods "Collation" subsection alongside the other source-database versions.
- The heuristic overlay rule set needs codifying as a repo script (probably a YAML or TSV file listing patterns and their direction assignments) before the paper's claim that overlays are "applied globally" is defensible. The electron-carrier universal per-electron anchor prototype (Fdx / Trx / Flx / Grx, see `Scripts/Thermodynamics/Update_Compound_GroupContribution_Energies.py` future commits) is one candidate rule-set component.
- Whether OpenTECR gets its own Methods paragraph or is folded into the eQuilibrator refresh is an open decision (see `PAPER_2026_PLAN.md` open decision #2).
- Quinone/quinol dGpredictor error audit needs execution before we can defensibly report that source's coverage numbers for reactions using these carriers.
