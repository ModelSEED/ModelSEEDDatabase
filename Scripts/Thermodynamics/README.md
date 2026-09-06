# ModelSEED Biochemistry Database Scripts (Thermodynamics)

Here we have the scripts used to handle the thermodynamics data that
came from both the Group Contribution approach (<a
href="https://doi.org/10.1529/biophysj.107.124784">Jankowski et
al. 2008</a>) and from eQuilibrator (<a
href="https://doi.org/10.1371/journal.pcbi.1003098">Noor et
al. 2013</a>). Our approach in handling this data is described in the
<a
href="https://www.biorxiv.org/content/10.1101/2020.03.31.018663v2">paper</a>.

## Order of execution

The general order is that the energies from the application of the Group Contribution (GC) approach
are stored in the database first, and then the energies from eQuilibrator (EQ), which, in most
cases, take precedence, are used to overwrite the energies in the database

### Additive per-method thermodynamics

Each reaction and compound keeps every method's estimate **additively** in its
`thermodynamics` dict rather than collapsing them into a single value:

```json
"thermodynamics": {
    "Group contribution": [4.15, 1.22, "="],
    "eQuilibrator":       [-3.46, 0.05, ">"],
    "dGPredictor":        [-3.82, 0.02, ">", 0.94]
}
```

Record shape differs by kind and by source:

| kind | source | record |
|---|---|---|
| reaction | Group contribution, eQuilibrator | `[dg, err, operator]` |
| reaction | dGPredictor | `[dg, err, operator, coverage]` |
| compound | Group contribution, eQuilibrator | `[dg, err]` |
| compound | dGPredictor | `[dg, err, coverage]` |

kcal/mol throughout. Compounds carry no operator — a formation energy has no
direction. dGPredictor's trailing `coverage` is appended, so positional readers
of `[0]`/`[1]`/`[2]` are unaffected; anything **rewriting** these lists must
preserve elements past the operator.

The **operator** (`>`, `<`, `=`, `?`) is that estimate's own thermodynamic
direction, computed from that source's own energy with **that source's own rule
set** — not the canonical reversibility. Passing `source=` to
`reversibility_from_energy` selects the set; omitting it silently falls back to
`GC`, scoring an eQuilibrator or dGPredictor energy with Group-Contribution
rules. The registry is documented in the next section.

`Add_Reaction_Thermodynamics_Operators.py` refreshes every operator from the
stored energies without needing the upstream inputs. Because each updater
already computes its operator at write time, a clean run reports
**"Entries refreshed/added: 0"**; a non-zero count means a record was written
with the wrong rule set.

These per-source records sit **next to** the canonical top-level
`deltag` / `deltagerr` / `reversibility` fields and never replace them.

> **The canonical fields are being retired** in favour of this dict. They are
> currently stale with respect to every source beneath them: derived before the
> 2026-08 eQuilibrator regeneration and dGPredictor retrain, and some promoted
> from `dGPredictor-ModelSEED`, which no longer exists.
> `Promote_Reaction_Thermodynamics_to_Canonical.py` never overwrites an existing
> canonical value, so re-running it does not repair them.

### Source-specific reversibility heuristics

The direction cascade lives in `reversibility_heuristics.py` and is selected
**per thermodynamic data source**, because the sources do not fail the same way.
`Estimate_Reaction_Reversibility.py` picks the rule set from its argument:

```
./Estimate_Reaction_Reversibility.py            # top-level deltag, GC rules
./Estimate_Reaction_Reversibility.py GC         # Group contribution, GC rules
./Estimate_Reaction_Reversibility.py EQ         # eQuilibrator energies, EQ rules
./Estimate_Reaction_Reversibility.py EQ --heuristics EQ2   # eQuilibrator 2.0 rules
./Estimate_Reaction_Reversibility.py EQ --heuristics GC    # pre-split behaviour
```

**`GC` is the default.** Any source without a rule set of its own — including
both dGPredictor variants — falls back to it, and the GC cascade is unchanged,
so `GC` and unfiltered runs still reproduce the historical report exactly.

| set | cascade | citation |
|-----|---------|----------|
| `GC` | ATPS → ABC transporter → MdeltaG bounds → mMdeltaG band → low-energy points → `=` | Jankowski et al. 2008 |
| `EQ` | ATPS → ABC transporter → undecomposable → uncorrected transport → ln Γ ± 1σ → `=` | Beber et al. 2022 over Noor et al. 2012 |
| `EQ2` | same, with the confidence margin off (ln Γ point estimate) | Flamholz et al. 2012 + Noor et al. 2012 |

The same registry drives the per-method operators, so
`thermodynamics['eQuilibrator'][2]` is now computed with the EQ rules while the
other methods keep the GC rules.

#### What the EQ rules do

Direction comes from the **reversibility index** Γ of Noor et al. 2012 — the
fold change every reactant concentration must undergo to reverse the reaction:

```
ln Γ = (2 / Σ|ν|) · ΔG′m / RT        ΔG′m = ΔG′° + RT · Σν · ln(10⁻³)
```

with water and protons excluded from both sums, matching
`equilibrator_cache.reaction.items(protons=False, water=False)`. A reaction is
irreversible when `|ln Γ| − σ > ln(1000)`, i.e. Noor's headline 3 µM–3 mM window
around 100 µM. `Scripts/Tests/test_eq_heuristics.py` verifies this reproduces
eQuilibrator's own published `ln_reversibility_index` (column 4 of
`MetaNetX_Reaction_Energies.tbl`) on 17,771 reactions.

Beber 2022 defines no reversibility index; what eQuilibrator 3.0 contributes is
rigorous uncertainty, and that is what the two gates ahead of the index rule
encode. Both address defects in the stored eQuilibrator data:

* **`EQ:undecomposable`** — a σ of ~1e5 kJ/mol (`RMSE_inf`) marks a reaction
  eQuilibrator could not decompose. For the majority of these it does not return
  a prediction at all: `GibbsEnergyPredictor.standard_dg` short-circuits with

  ```python
  if residual:
      return Q_(0, "kJ/mol").plus_minus(self.preprocess.RMSE_inf)
  ```

  i.e. **literally zero**, and the nonzero ΔG′° we store is only the Legendre/pH
  transform that `standard_dg_prime` adds on top of that zero. Of the 4,607
  σ-flagged reactions in the table, **2,895** have σ exactly `RMSE_inf` and are
  this zeroed case — 61 of them are still exactly `0.000000` because they have no
  net proton change for the transform to act on. The remaining **1,712** carry
  `RMSE_inf * ||sigma_inf||` (e.g. 1/√2) and do have a real mean, but with an
  undetermined direction in the covariance.

  Either way the value cannot support a directional call, so the gate returns
  `?`. Previously the GC bounds rule could never fire at that width and these
  fell through to a permissive `=`. (Real σ tops out at 65.35 kcal/mol against a
  marker of 23,900.57, so the cut at 1e4 kJ/mol sits in an empty gap.)
* **`EQ:transport-uncorrected`** — `Retrieve_eQuilibrator_Reactions_Energies.py`
  keys its MetaNetX formula on compound id and therefore discards compartment,
  so any species on both sides nets out; 1,102 transport reactions carry a ΔG′°
  for a *different* reaction. Independently, Beber 2022 notes the transformed
  framework needs a `−N_H·RT·ln(10^ΔpH) − Q·FΔΦ` term across a membrane that
  this pipeline never applies. ATP synthase and ABC transporters are still
  decided structurally, ahead of this gate.

Reactions whose ln Γ interval merely straddles the threshold are reported as
`EQ:ambiguous` and called `=`: Noor treats Γ as a continuous index and reserves
the directional call for clear cases.

#### Energy source for the EQ run

`Estimate_Reaction_Reversibility.py EQ` reads
`thermodynamics['eQuilibrator']`'s own dG and σ, **not** the top-level `deltag`.

Reading `deltag` was correct as originally designed. From 2019 (`c263e233`)
through 2023-09-11 (`17c9739b`), `Update_Reaction_eQuilibrator_Energies.py`
overwrote the canonical energy and stamped the note that gates this step:

```python
reactions_dict[rxn]['deltag']    = float(eq_reactions[rxn]['dg'])
reactions_dict[rxn]['deltagerr'] = float(eq_reactions[rxn]['dge'])
if('EQU' not in notes_list): notes_list.append('EQU')
```

so `deltag` *was* the eQuilibrator value wherever `EQU` was set, and
`DB_LEVEL_NOTE = {"EQ": "EQU"}` gated it exactly right.

**`3e50646b` (2023-09-13) removed the write half and left the read half in
place.** The updater moved to `thermodynamics['eQuilibrator']` only, stopping
both the `deltag` overwrite and the `EQU` stamp — but this step kept reading
`deltag`. Later GC rebuilds (Convention A, `ad34d6ab`) and the PR #265
promotion then rewrote `deltag` from other sources. Today only **1,797 of
25,028** reactions with an eQuilibrator record have
`deltag == thermodynamics['eQuilibrator'][0]`; restricted to the 17,094 still
carrying `EQU`, only **1,586** do. That note is now a fossil of the pre-2023
pipeline. (`EQU` is not `EQC`: the retrieval script writes `EQC`/`EQP` to mean
"all reagents had eQuilibrator structures", a different claim.)

Pointing the EQ run at the eQuilibrator sublist restores the original intent.
The EQ rules also need eQuilibrator's own σ, which the top-level `deltagerr`
never carries.

Consequence: canonical `reversibility` for these reactions now derives from the
eQuilibrator energy while canonical `deltag` mostly remains the GC value.
Reconciling `deltag` itself is the job of
`Promote_Reaction_Thermodynamics_to_Canonical.py`, which is deliberately not
part of `Rerun_Thermodynamics.sh`.

### The reversibility index as a standalone callable

`reversibility_index.py` is the Noor et al. 2012 / eQuilibrator 2.0 index
(Bioinformatics 28:2037; NAR 40:D770) on its own, with no cascade around it:

```python
from reversibility_index import ln_reversibility_index, direction_from_index

ln_gamma = ln_reversibility_index(rxn["stoichiometry"], dg)   # kcal/mol in
operator = direction_from_index(ln_gamma)                     # '>', '<' or '='
```

Γ is the fold change every reactant concentration would have to undergo,
symmetrically, before the reaction ran backwards; `ln Γ = (2/Σ|ν|)·ΔG′m/RT`,
with water and protons dropped from both coefficient sums and every reagent
placed at 1 mM. `|ln Γ| > ln(1000)` is Noor's cut.

`Context.ln_gamma` delegates here, so the cascade and any standalone caller
cannot drift apart. Validated against eQuilibrator's own published
`ln_reversibility_index` column in `MetaNetX_Reaction_Energies.tbl`: **17,776
reactions agree**, and 93% of the residual mismatches are the compartment-
collapse defect described below.

**Transport reactions get no special treatment anywhere in this module.** Every
`(compound, compartment)` pair is its own species; there is no membrane term, no
compartment collapse and no structural shortcut. That is what lets a transport
call be compared with a cytosolic one on the same axis.

#### What the index says about ATP synthase and ABC transporters

`Scripts/Tests/test_reversibility_index.py --report` scores both families as
ordinary chemistry. Two results worth knowing before trusting the index on them:

* **ATP synthase (15 reactions).** All 15 come out with `|ln Γ| = 11.97` —
  *identical*, because the translocated protons cancelled in the collapsed
  MetaNetX formula, so what eQuilibrator actually scored was
  `ADP + Pi ⇌ ATP + H₂O` every time. The index therefore calls the most
  famously reversible enzyme in the cell irreversible, and splits 7 `>` / 8 `<`
  purely on how each reaction happened to be written. Letting the protons count
  (`exclude_protons=False`) returns 13 of 15 to `=`.
* **ABC transporters (1,150 scorable).** 94.1% agreement with the structural
  ATP-sign rule, 29 outright flips, 39 called reversible. The agreement is
  thinner than it looks: the transported substrate also cancels in the collapsed
  formula, leaving plain ATP hydrolysis at −6.54 kcal/mol and `ln Γ = −7.18`
  against a cut of −6.91. 1,074 of the 1,111 hard calls sit within 15% of the
  threshold.

Both families are decided by structural shortcuts in the `GC` and `EQ` sets
precisely because the energy behind them cannot be trusted. The index does not
change that; it just makes the untrustworthiness measurable.

### Applying the index to dGPredictor

`SOURCE_HEURISTIC_SET` now routes `dGPredictor` and `dGPredictor-ModelSEED` to
`DGP_HEURISTICS` (`make_ri_heuristics(z=1.0, sigma_gate=None)`) instead of
letting them fall through to the GC concentration bounds. Neither source has a
"could not decompose" marker, so there is no sentinel to gate on; the one-sigma
margin is kept because silent extrapolation, not a loud refusal, is how these
models fail.

The margin lands very differently on the two, and that is the point:

| source | σ median | σ max | `=` from the margin | vs GC operators |
|---|---:|---:|---:|---|
| `dGPredictor` | 0.35 | 6.10 | 784 / 27,715 | 2,816 re-scored |
| `dGPredictor-ModelSEED` | 21.17 | 2,039.14 | 19,512 / 31,924 | — |

Use the `RI` rule set (z = 0, the index exactly as published) to see the
point-estimate answer for either source.

### Comparing rule sets on one set of energies

`Compare_Reversibility_Heuristics.py` hands two or more rule sets the *same*
`(dG, σ)` pair, so disagreement is attributable to the rules alone:

```bash
./Compare_Reversibility_Heuristics.py --sets GC EQ EQ2 RI
./Compare_Reversibility_Heuristics.py --source dGPredictor --sets GC DGP
```

On eQuilibrator energies, GC and EQ **agree on 71.6%** of 25,028 reactions. The
disagreement is dominated by one cell: **4,369 reactions GC calls `=` that EQ
calls `?`**, which is the undecomposable population the GC bounds rule swallows.
The next largest is 1,447 that GC calls `>` and EQ calls `=`. Notably, on GC's
side the terminal `default` rule fires 8,556 times and the concentration-bounds
rules 8,714 — so a third of GC's answers on this source come from a rule that
never looks at σ.

### Checking eQuilibrator's error returns

`Check_eQuilibrator_Energy_Errors.py` looks for the records where eQuilibrator
reported a failure as a number. It never raises, so these are otherwise silent:

```bash
./Check_eQuilibrator_Energy_Errors.py --table --tsv findings.tsv
./Check_eQuilibrator_Energy_Errors.py --self-test    # proves each check fires
```

| code | database | raw table | what it is |
|---|---:|---:|---|
| `UNDECOMPOSABLE` | 4,934 | 4,607 | σ at or past `RMSE_inf` = 1e5 kJ/mol — eQuilibrator declined |
| `ZERO_RETURN` | 632 | 184 | dG exactly 0.0 — the discarded-mean branch showing through |
| `COLLAPSED_FORMULA` | — | 1,178 | published `ln_RI` disagrees with the index recomputed from our stoichiometry |

`COLLAPSED_FORMULA` is database-invisible because column 4 of the table is
discarded on the way in (`_thermo_helpers.parse_two_col_energy_table` reads
columns 1 and 2 only), which is why `--table` exists.

`--self-test` runs 15 synthetic records through the same `check_record()` the
scans use — including negative controls that must stay clean — and asserts the
two headline cases explicitly: the 0-return energy and the 1e5 sentinel σ.

### Translating eQuilibrator's sentinel into ModelSEED's

`Normalize_eQuilibrator_Sentinels.py` rewrites eQuilibrator's failure returns in
ModelSEED's convention, so a consumer that correctly skips one skips the other.
The two are not compatible:

| | "no value" looks like |
|---|---|
| ModelSEED | `[10000000.0, 10000000.0, '?']` — unmistakable |
| eQuilibrator | a plausible signed energy with σ = 1e5 kJ/mol — indistinguishable |

```bash
./Normalize_eQuilibrator_Sentinels.py                 # dry run, report only
./Normalize_eQuilibrator_Sentinels.py --apply         # rewrite the JSON
./Normalize_eQuilibrator_Sentinels.py --table out.tbl # normalise at ingest instead
./Normalize_eQuilibrator_Sentinels.py --self-test
```

Detection is `check_record()` imported from `Check_eQuilibrator_Energy_Errors.py`,
so there is one definition of "sentinel" and not two. Nothing is written without
`--apply`.

**Scope: 4,934 of 25,028 records (19.7%)** — 4,590 `UNDECOMPOSABLE` plus 344 that
are also exactly zero. Only `thermodynamics['eQuilibrator']` is touched; the
sibling per-source triples and the flat `deltag`/`deltagerr` are left alone.
4 of the 4,934 are *also* the current canonical `deltag`; those are reported, not
rewritten, because choosing a replacement is
`Promote_Reaction_Thermodynamics_to_Canonical.py`'s job.

#### Why `zero_dg` is off by default

It is tempting to also blank the 288 records whose dG is exactly `0.0` with a
credible σ, reading them as the residual branch showing through where the
transform was nil. **They are not that branch** — the residual branch always
stamps `RMSE_inf` on σ, so by construction these are something else. 280 of them
carry `σ == 0.0` *exactly*, and 232 of those are one-in/one-out isomerizations:
L-lysine ⇌ D-lysine, D-threo- ⇌ D-erythro-isocitrate, 16α- ⇌ 16β-hydroxysteroid.
eQuilibrator's decomposition is stereo-blind, so both sides decompose to
identical groups and the difference is exactly zero with exactly zero propagated
error. That is a real (if uninformative) statement that the two stereoisomers
are isoenergetic, not a failure report, and blanking it throws it away.

Enable with `--classes undecomposable non_finite zero_dg`. If what you actually
want is the 8 records where a zero energy is genuinely suspicious, add
`--zero-dg-nonzero-sigma` to spare the stereo cases.

`COLLAPSED_FORMULA` is deliberately excluded: a wrong-but-finite energy from our
own compartment-dropping retrieval step is a different defect needing a
different fix (rebuild the formula), not a sentinel.

The underlying thermodynamics data is kept in
`../../Biochemistry/Thermodynamics`, one directory per source, each with its own
README recording provenance:

| directory | source | README |
|---|---|---|
| `ModelSEED/` | Group contribution — MFAToolkit MolAnalysis tables | — |
| `eQuilibrator/` | component contribution from ModelSEED structures | yes |
| `dGPredictor/` | retrained fragment model + staged predictions | yes |

`Retrieve_eQuilibrator_{Compound,Reactions}_Energies.py` built the superseded
MetaNetX-mediated tables and are kept for reference only; they no longer feed
the database. See `../../Biochemistry/Thermodynamics/eQuilibrator/README.md`.

If the underlying thermodynamics data hasn't changed, running the pipeline
should produce no diff:

```
./Rerun_Thermodynamics.sh
```

That script documents the full order, including the canonical-field steps that
are deliberately **not** run while those fields are being retired.
