# Source-specific reaction-reversibility heuristics

**Branch:** `eQuilibrator-fix` (from `ModelSEED/ModelSEEDDatabase` `dev` @ `49563c6f`)
**PR:** [ModelSEED/ModelSEEDDatabase#285](https://github.com/ModelSEED/ModelSEEDDatabase/pull/285) — `Cooper-Taylor:eQuilibrator-fix` → `ModelSEED:dev`
**Date:** 2026-08-18

---

## 1. The question we started from

ModelSEED assigns every reaction a direction (`>`, `<`, `=`, `?`) from its stored
free energies. We wanted that assignment to depend on *which* thermodynamic method
supplied the energy, rather than applying one rule set to all of them — and
specifically to give eQuilibrator a cascade built on eQuilibrator's own methodology.

## 2. What we found in `dev` first

Before writing anything we established what `dev` actually does. The answer:
**the Group-Contribution heuristics decide everything.** The eQuilibrator
reversibility index is never used, as a cutoff or otherwise.

Evidence, all verified against `origin/dev` rather than inferred:

- `DEFAULT_HEURISTICS` is the only cascade defined, and it is passed literally at
  both call sites in `Estimate_Reaction_Reversibility.py`.
- `make_ln_reversibility_index_heuristic` exists at `reversibility_heuristics.py:264`
  but is **never invoked** — every apparent use is either the module docstring
  showing how one *could* wire it up, or a back-compat re-export. `LN_RI_THRESHOLD`
  sits in the file purely as that function's unused default argument.
- `Retrieve_eQuilibrator_Reactions_Energies.py` calls `ln_reversibility_index()` and
  writes it as **column 4** of `MetaNetX_Reaction_Energies.tbl` for 23,592 reactions.
  `parse_two_col_energy_table` reads only columns 1–2. Column 4 is computed, stored,
  and discarded on every run.
- Executing `dev`'s own cascade against `dev`'s own JSONs confirms it: every deciding
  rule in both the `GC` and `EQ` runs comes from the GC vocabulary
  (`MdeltaG(Max)`, `mMdeltaG`, `lowE`, `MdeltaG(Min)`, `ABCT`, `ATPS`, `default`).

### Three defects this surfaced

**(a) The `EQ` run never read eQuilibrator energies.** `_energy_for(entry, 'EQ')`
pulls the canonical `deltag`, merely *gated* on eQuilibrator eligibility. Since the
additive-thermodynamics refactor, no updater overwrites `deltag` — so only **1,797 of
25,028** reactions with an eQuilibrator record actually had
`deltag == thermodynamics['eQuilibrator'][0]`. The other **23,140** were scored on the
Group-Contribution number and then labelled eQuilibrator. Canonical `deltagerr` is
*never* the eQuilibrator sentinel (0 cases), which confirms it is the GC error
throughout.

**(b) eQuilibrator's "cannot decompose" marker was read as an error bar.** A σ of
~1e5 kJ/mol is a flag; **4,933** reaction records carry it. The GC bounds rule cannot
fire at that width, so those reactions fell through to the narrow mM band or the bare
`default` and came out `=`. Observed *real* σ tops out at **65.35** kcal/mol against a
marker of **23,900.57**, so the two populations are cleanly separated.

**(c) Transport energies are for the wrong reaction.**
`Retrieve_eQuilibrator_Reactions_Energies.py` builds its formula keyed on MetaNetX id,
discarding compartment, so any species on both sides nets out. **1,102** transport
reactions carry a ΔG′° for a different reaction — `rxn12518` comes out at
+524 kcal/mol. A further **76** non-transport reactions are hit by the same collapse
via stereo-neutral InChIKey matching (e.g. `rxn00816`, where D-glucose and galactose
merge into one compound).

## 3. What eQuilibrator actually provides

Worth recording, because it constrains what a faithful EQ cascade can be.

**eQuilibrator ships no directionality classifier.** `ComponentContribution` exposes
one directionality-adjacent primitive — `ln_reversibility_index` — and it returns a
number, not a verdict. `PhasedReaction` has nothing; "reversible" appears there only as
a direction-insensitive hashing flag. Beber 2022 defines no reversibility index at all
and explicitly hands directionality constraints to downstream tools (multiTFA, pyTFA,
PTA), offering MDF and ECM instead — both of which *score* pathways rather than predict
direction.

The index itself is Noor et al. 2012, used by eQuilibrator 2.0 (Flamholz 2012):

```
ln Γ = (2 / Σ|ν|) · ΔG′m / RT          ΔG′m = ΔG′° + RT · Σν · ln(10⁻³)
```

with water and protons excluded from both sums. Noor deliberately declines to give a
cutoff — Γ is presented as a *continuous* index precisely to replace binary annotation
— but states a headline window of γ = 1000 (3 µM–3 mM around 100 µM), under which
~55% of reactions are reversible.

What Beber 2022 (eQuilibrator 3.0) contributes: ~10,000 → ~1,000,000 compounds, a move
from KEGG to MetaNetX identifiers (which is why ModelSEED ids map at all), component
contribution throughout, and — the headline — **full covariance** uncertainty via
`standard_dg_prime_multi`, on the argument that treating reaction uncertainties as
independent "usually overestimates the uncertainty and may violate the first law."

**So `EQ_HEURISTICS` is not a ported eQuilibrator cascade — no such thing exists.**
Its provenance is mixed, and honestly so:

| component | provenance |
|---|---|
| ln Γ formula, exclusions, 1 mM convention | eQuilibrator's own implementation (verified) |
| ln(1000) threshold | Noor 2012's stated window — their convention, not their cutoff |
| ±1σ confidence margin | **ours**, motivated by Beber 2022 |
| undecomposable + transport gates | **ours**, from eQuilibrator's documented failure modes |

## 4. What we built

A rule-set registry in `reversibility_heuristics.py`, selected per source:

| set | cascade | assigned to |
|---|---|---|
| `GC` **(default)** | ATPS → ABCT → MdeltaG bounds → mMdeltaG band → lowE → `=` | Group contribution, dGPredictor, dGPredictor-ModelSEED, anything unrecognised |
| `EQ` | ATPS → ABCT → undecomposable → uncorrected transport → ln Γ ± 1σ → `=` | eQuilibrator |
| `EQ2` | same, margin off (bare point estimate) | opt-in via `--heuristics EQ2` |

Structural rules (ATP synthase, ABC transporter) lead in both sets because they need no
energy, so neither eQuilibrator data defect can reach them.

Ambiguous reactions — where the ln Γ interval straddles the threshold — are labelled
`EQ:ambiguous` and called `=`, matching the GC cascade's permissive terminal rule and
Noor's reservation of directional calls for clear cases.

**The GC cascade is byte-for-byte unchanged.** `DEFAULT_HEURISTICS` remains as an alias
so existing importers keep working.

Other changes: `Estimate_Reaction_Reversibility.py` gained a `--heuristics` override
and now sources EQ energy from the eQuilibrator sublist; `_thermo_helpers`,
`Add_Reaction_Thermodynamics_Operators` and `Promote_*` now pass the source label
through, so each stored operator uses its own rule set.

## 5. Results

Per-method operators, vs `dev`:

| source | records | `>` | `<` | `=` | `?` | changed |
|---|---:|---:|---:|---:|---:|---:|
| Group contribution | 56,002 | 9,579 | 1,413 | 16,321 | 28,689 | **0** |
| eQuilibrator | 25,028 | 8,655 | 1,440 | 9,791 | 5,142 | **7,103** |
| dGPredictor | 27,715 | 13,125 | 3,134 | 11,456 | 0 | **0** |
| dGPredictor-ModelSEED | 31,924 | 9,922 | 1,231 | 20,771 | 0 | **0** |

Only eQuilibrator moves. The three GC-scored columns being byte-identical is the
control: it isolates the change to the new cutoff rather than the refactor.

Canonical reversibility:

| | dev | now | Δ |
|---|---:|---:|---:|
| `>` | 12,670 | 9,584 | −3,086 |
| `<` | 2,323 | 1,603 | −720 |
| `=` | 15,059 | 14,002 | −1,057 |
| `?` | 25,960 | 30,823 | **+4,863** |

Deciding rule in the EQ run: `Incomplete` 24,828 · **`EQ:lnGamma` 8,944** ·
`EQ:reversible` 8,905 · `Incomplete (GCC)` 6,146 · `EQ:undecomposable` 4,933 ·
`ABCT` 1,151 · `EQ:ambiguous` 871 · `EQ:transport-uncorrected` 209 · `ATPS` 15 ·
`Empty` 10.

### Reading it

**The database gets less confident, and that is the finding.** The +4,863 unknowns
track the 4,933 undecomposable records almost exactly. Those reactions had been
receiving a permissive `=` or an inherited direction on the strength of a ΔG′° carrying
±23,900 kcal/mol — an unread sentinel, not evidence.

The bare no-evidence `default` fallback went from deciding **7,120** reactions to
**zero**: everything that used to land there now gets a real index call, an explicit
`EQ:reversible`, or an honest `?`.

The **3,086 lost `>` calls** are what to scrutinise before merging. Two populations are
mixed in there: reactions where the index genuinely disagrees with the GC bounds rule
on real numbers, and reactions whose only evidence was the sentinel. They are separable
from the generated report.

## 6. Verification

- **New** `Scripts/Tests/test_eq_heuristics.py` — the ln Γ implementation reproduces
  eQuilibrator's own published `ln_reversibility_index` on **17,771** reactions. Every
  residual is a MetaNetX-collapsed reaction, asserted rather than assumed. Plus unit
  assertions per rule and on the registry defaulting to GC.
- `Scripts/Tests/test_reaction_direction.py` — now separates invariant sources (GC,
  dGPredictor: must match exactly) from intentionally re-scored ones (eQuilibrator),
  with `--strict` to require equality everywhere. Passes.
- `Rerun_Thermodynamics.sh` is **idempotent** — a second run reproduced all 61 JSONs
  byte-for-byte (`Entries refreshed/added: 0`).
- Spot checks: `rxn00001` → `>` (ln Γ −9.18, pyrophosphatase) · `rxn00017` → `?`
  (sentinel σ) · `rxn12518` → `?` (collapsed transport) · `rxn08173` ATP synthase → `=`
  structurally.

## 7. Deliberately out of scope

- **Canonical `deltag` is untouched** and mostly remains the Group-Contribution value,
  so `reversibility` and `deltag` can now disagree in provenance. Reconciling them is
  `Promote_Reaction_Thermodynamics_to_Canonical.py`, which is not part of
  `Rerun_Thermodynamics.sh` and was not run.
- **No dGPredictor estimate reaches canonical `reversibility`.**
  `Estimate_Reaction_Reversibility.py DGP` works, but the rerun script only invokes GC
  and EQ, so 27,715 + 31,924 operators sit in `thermodynamics` influencing nothing.
  With `?` now at 30,823 that is a real coverage opportunity, but promoting it is a
  policy call.
- **Travis is already inert** — `.travis.yml` validates `Biochemistry/reactions.json`
  and `compounds.json`, neither of which exists or is tracked. Not made worse here.

## 8. Suggested follow-ups

1. Replace the σ-threshold gate with **`is_using_group_contribution()`** — eQuilibrator's
   own verdict on whether it could decompose a reaction, removing our magic number.
2. Use **`standard_dg_prime_multi()`** for covariance. We currently store 25,028
   independent scalar σ's, which Beber 2022 says systematically *overestimates*
   uncertainty — directly inflating the ±1σ margin in the EQ cutoff.
3. Use **`multicompartmental_standard_dg_prime()`** for the 1,102 transport reactions
   currently gated to `?`. Needs real ΔΦ and ΔpH per compartment pair, but converts a
   discard into an estimate.
4. Fix the retrieval script's MetaNetX collapse at source, so transport and
   stereoisomer reactions get the reaction they actually are.
