# How eQuilibrator's undecomposable reactions became confident directions

**Forensic walkthrough of the defect fixed by [PR #285](https://github.com/ModelSEED/ModelSEEDDatabase/pull/285).**
Code references are to `origin/dev` @ `49563c6f` unless stated otherwise.
Companion documents:
`equilibrator_reversibility_heuristics_2026-08-18.md` (what was built and why) and
`equilibrator_deltag_read_write_defect_2026-08-18.md` (the separate `deltag` read/write
defect, which also dates to `3e50646b`).

---

## Summary

For **4,934 reactions**, eQuilibrator explicitly declines to estimate a free energy. It
signals this by returning **zero** with an uncertainty of 1e5 kJ/mol. A Legendre/pH
transform is then added on top of that zero — by eQuilibrator, inside the same API call
— so what arrives at our pipeline is a plausible-looking nonzero number carrying an
enormous error bar.

On `dev` those 4,934 reactions receive **confident directions**: 2,940 get a hard `>` or
`<`, 1,903 get `=`, and only 91 get `?`. The directional calls come from a rule that
reads the *sign* of the energy — which for these reactions is the sign of the transform.
**Reactions are being called irreversible on the basis of how many protons they release
at pH 7.**

The behaviour was introduced on **2023-09-13** by commit `3e50646b` and has been live
ever since. A **2026-07-05** eQuilibrator rerun expanded coverage into a pipeline that no
longer had the guard, turning a handful of such records into thousands.

---

## 1. The pipeline

### Step 1 — Retrieve

`Scripts/Thermodynamics/Retrieve_eQuilibrator_Reactions_Energies.py` maps ModelSEED
compounds to MetaNetX ids, builds a reaction formula, and makes **one call**:

```python
result = equilibrator_calculator.standard_dg_prime(equilibrator_reaction)
dG0_prime   = str(result.value.to('kilocal / mole').magnitude)   # kJ -> kcal only
uncertainty = str(result.error.to('kilocal / mole').magnitude)
ln_RI       = equilibrator_calculator.ln_reversibility_index(equilibrator_reaction)
output_handle.write("\t".join([rxn, dG0_prime, uncertainty, ln_RI]) + "\n")
```

Conditions are set at pH 7.0, I = 0.25 M, T = 298.15 K. Our only arithmetic is the
kJ→kcal unit conversion.

### Step 2 — What eQuilibrator does inside that call

Two methods in `component_contribution/predict.py`:

```python
def standard_dg_prime(self, reaction, p_h, ionic_strength, temperature, p_mg):
    return self.standard_dg(reaction) + reaction.transform(       # transform ALWAYS added
        p_h=p_h, ionic_strength=ionic_strength,
        temperature=temperature, p_mg=p_mg,
    )

def standard_dg(self, reaction):
    mu, sigma_fin, sigma_inf, residual = self.get_reaction_prediction(reaction)
    if residual:
        return Q_(0, "kJ/mol").plus_minus(self.preprocess.RMSE_inf)   # mean DISCARDED
    else:
        std = np.linalg.norm(sigma_fin, 2) \
            + self.preprocess.RMSE_inf * np.linalg.norm(sigma_inf, 2)
        return Q_(mu, "kJ/mol").plus_minus(max(self.preprocess.RMSE_eps, std))
```

When a reaction has a component outside both the reactant- and group-contribution spans,
`standard_dg` returns **literally zero** — not a minimum-norm solution, not a projection.
The computed mean is thrown away. `standard_dg_prime` has **no residual check of its
own**; it adds the transform unconditionally. So we receive `0 + transform`.

The transform is protonation and Mg²⁺ bookkeeping (Θ_H, Θ_Mg summed over reagents at the
given conditions). It is computable for any compound with known pKa's and has **nothing
to do with whether the reaction's chemistry is predictable**. That is why the output
looks like a real energy.

Two populations result, both flagged by σ:

| population | in table | in JSON | σ | what the number is |
|---|---:|---:|---|---|
| residual branch | 2,895 | 3,037 | exactly `RMSE_inf` = 1e5 kJ/mol | **zero + transform** |
| undetermined direction | 1,712 | 1,897 | `RMSE_inf · ‖σ_inf‖` (e.g. 1/√2) | real mean, inflated error |

**Direct confirmation.** Among residual-branch reactions with *no net proton change* —
nothing for the transform to act on — 35% come out below 0.01 kcal/mol and **61 are
exactly `0.000000`**:

```
rxn02677  dG'o = 0.000000000 kcal/mol   (1) Cycloeucalenol[0] <=> (1) Obtusifoliol[0]
rxn04059  dG'o = 0.000000000 kcal/mol   (1) 16-Epivellosimine[0] <=> (1) Vellosimine[0]
rxn04757  dG'o = 0.000000000 kcal/mol   (1) Violaxanthin[0] => (1) Neoxanthin[0]
```

The zero shows through unmasked. `rxn00017` releases 2 H⁺, so its transform is
−20.46 kcal/mol — and that is the *entire content* of the number we store for it.

Real σ across the whole table tops out at **65.35** kcal/mol against a marker of
**23,900.57**, so the two regimes are cleanly separated with nothing in between.

### Step 3 — Store

`Update_Reaction_eQuilibrator_Energies.py` on `dev` is a three-line delegation to
`_thermo_helpers`:

* `_thermo_helpers.py:120-136` — `parse_two_col_energy_table` reads `array[1]` and
  `array[2]` only. **Column 4, the reversibility index, is read into `array[3]` and
  discarded on every run.**
* `_thermo_helpers.py:413` — `run_reaction_lookup_update` calls
  `_per_source_operator(reactions_dict[rxn], dg, dge)` with **no label argument**, so
  every source is scored with Group-Contribution rules.

There is **no error threshold anywhere** in this path. The σ itself is stored faithfully;
nothing is lost at this stage.

### Step 4 — Decide direction

`Estimate_Reaction_Reversibility.py:75`, inside `estimate_one`:

```python
status, thermoreversibility, source_label = run_reversibility(
    rxn_entry, top_level_energy(db_level), DEFAULT_HEURISTICS)
```

Two things are hardcoded here:

1. **`top_level_energy`** resolves via `reversibility_heuristics.py:80` —
   `rxn_dg = rxn_entry['deltag']`, the *canonical* field, not the eQuilibrator sublist.
2. **`DEFAULT_HEURISTICS`** (`reversibility_heuristics.py:254`) is the GC cascade,
   regardless of `db_level`.

Against σ = 23,900, that cascade behaves as follows:

| rule | line | behaviour |
|---|---|---|
| `stored_bounds_heuristic` | 215 | reads `dge`; `hi > 0` and `lo < 0` always → returns `None` |
| `mmdeltag_band_heuristic` | 235 | **σ not referenced** |
| `low_energy_heuristic` | 241 | **σ not referenced** — takes the *sign* of `mMdeltaG` |
| `default_heuristic` | 249 | always `=` |

σ *is* read — once. An enormous value degrades that rule to "no opinion", and control
passes to three rules that cannot form one. The result is written to canonical
`reversibility` at `Estimate_Reaction_Reversibility.py:163`.

---

## 2. When this started: `3e50646b`

**Authored and committed 2023-09-13**, "Updating integration of thermodynamics data in
biochemistry for reactions". The preceding commit, `17c9739b`, was two days earlier.

### The last functional version — `17c9739b` (2023-09-11)

```python
    #Here we establish an arbitrary threshold of 100 for the error, if the error
    #is too big, we don't use it

    if(float(eq_reactions[rxn]['dge']) > 100):
        continue                                                      # (A) guard

    notes_list=reactions_dict[rxn]['notes']
    if(not isinstance(notes_list,list)):
        notes_list=list()

    reactions_dict[rxn]['deltag']=float(eq_reactions[rxn]['dg'])       # (B)
    reactions_dict[rxn]['deltagerr']=float(eq_reactions[rxn]['dge'])   # (B)
    if('EQU' not in notes_list):
        notes_list.append('EQU')                                       # (C)
    reactions_dict[rxn]['notes']=notes_list
```

### The first broken version — `3e50646b` (2023-09-13)

```python
    #Here we establish an arbitrary threshold of 100 for the error, if the error
    #is too big, we shouldn't use it, but for storing the database we keep it now
    if(float(eq_reactions[rxn][1]) > 100):
        pass                                                          # (A) neutered

    # values always saved as list of energy and error
    if(not isinstance(reactions_dict[rxn]['thermodynamics'],dict)):
        reactions_dict[rxn]['thermodynamics'] = dict()
    if(label not in reactions_dict[rxn]['thermodynamics']):
        reactions_dict[rxn]['thermodynamics'][label]=list()
    reactions_dict[rxn]['thermodynamics'][label]=eq_reactions[rxn]
                                                                      # (B),(C) deleted
```

### Three separate breakages in one commit

**(B) The `deltag` / `deltagerr` overwrite was deleted.** This is the one that matters
most. `Estimate_Reaction_Reversibility.py` still reads `rxn_entry['deltag']`. The write
half was removed; the read half was not updated. **The original design was correct** —
`deltag` genuinely was the eQuilibrator value wherever the note was set, so reading it
was reading eQuilibrator energies exactly as intended. Only the write side broke.

**(C) The `EQU` stamp was deleted.** `DB_LEVEL_NOTE = {"EQ": "EQU"}` still gates on it.
`EQU` means "this deltag came from eQuilibrator", written by the *updater* — distinct
from `EQC`, which the *retrieval* script writes to mean "all reagents had eQuilibrator
structures". Both survive in the data (27,024 `EQC`, 17,094 `EQU`) but `EQU` stopped
being maintained, so it is now a 2023 fossil: only **1,586** of those 17,094 still have
`deltag` equal to the eQuilibrator value.

**(A) `continue` became `pass`.** `pass` is a no-op, so the branch has no effect. The
comment change states the intent openly — *"we shouldn't use it, but for storing the
database we keep it now"* — so storing everything was deliberate. The defect is that the
threshold was never reinstated at any point of *consumption*. Under `17c9739b`,
`continue` meant an undecomposable reaction was skipped entirely and kept its
Group-Contribution direction. Afterwards every record flows into
`thermodynamics['eQuilibrator']` and `Add_Reaction_Thermodynamics_Operators.py` computes
an operator for all of them.

`Promote_Reaction_Thermodynamics_to_Canonical.py` later reintroduced `MAX_ERR = 100.0`,
so *promotion* is protected. The reversibility cascade never got the guard back.

### It has been live ever since

The file was superseded, but **every later version preserved the behaviour**:

| commit | date | writes `deltag`? | stamps `EQU`? | error guard |
|---|---|---|---|---|
| `17c9739b` | 2023-09-11 | yes | yes | `continue` |
| **`3e50646b`** | **2023-09-13** | **no** | **no** | `pass` (dead) |
| `c887ede8` | 2026-06-04 | no | no | none |
| `ad34d6ab` | 2026-08-07 | no | no | none |

`c887ede8` ("Add additive per-source thermodynamic descriptions") refactored *around* the
defect rather than fixing it; `ad34d6ab` left it alone. The behaviour has therefore been
continuously live on `dev` from **2023-09-13 to 2026-08-18** — just under three years —
and is still present at `49563c6f`.

It reached `dev` through the merges of **2024-01-16**, **2024-05-02** and **2024-12-11**,
and has been the shipped behaviour for every downstream consumer since.

### Why 2026 made it far worse

`MetaNetX_Reaction_Energies.tbl` regeneration history:

```
32cd10c4  2026-07-05  Re-run eQuilibrator on ModelSEED structures (fix API drift; expand coverage)
6336bd44  2020-04-24  Removing energies for reactions with SRUs
a0890dea  2019-09-23  Updating reaction energies and notes following updates
f3f652f9  2019-07-01  Adding reaction energies as retrieved from eQuilibrator using version 0.2.0
```

The prior regeneration was **2020-04-24**, when the `> 100 → continue` guard was still
live and undecomposable reactions were simply skipped. The **2026-07-05** rerun expanded
coverage substantially — and landed into a pipeline that no longer had the guard.

That is the causal chain: **2023 removes the guard and the `deltag` write; 2026 re-runs
eQuilibrator against a much larger compound set, generating far more undecomposable
reactions; nothing downstream filters them.**

---

## 3. Impact

What `dev` assigns to those 4,934 reactions, versus this branch:

| | dev | fixed |
|---|---:|---:|
| **eQuilibrator operator** `>` | 613 | 1 |
| `<` | 15 | 0 |
| `=` | 4,306 | 0 |
| `?` | **0** | 4,933 |
| **canonical reversibility** `>` | 2,262 | 1 |
| `<` | 678 | 0 |
| `=` | 1,903 | 0 |
| `?` | **91** | 4,933 |

Not one received `?` on `dev`. Reconstructing which GC rule produced each answer:

```
default   -> '='  3,642     the bare no-evidence fallback
mMdeltaG  -> '='    664
lowE      -> '>'    612     <- direction taken from a discarded zero
lowE      -> '<'     15
ABCT      -> '>'      1
```

The **627 directional calls from `lowE`** are the most damaging: that rule reads the sign
of `mMdeltaG`, which for the residual population is the sign of the Legendre transform.
In canonical `reversibility` the directional total is **2,940**, because the canonical
field is computed from `deltag` rather than the eQuilibrator sublist.

---

## 4. The fix

One rule, third in the eQuilibrator cascade:

```python
# Scripts/Thermodynamics/reversibility_heuristics.py:389
def eq_undecomposable_heuristic(ctx):
    if abs(ctx.dge) >= EQ_UNDECOMPOSABLE_SIGMA:      # 1.0e4 / 4.184 = 2390.06 kcal/mol
        return f"EQ:undecomposable: {ctx.dge:.0f}", "?"
    return None
```

The threshold sits in the empty gap between the largest genuine σ (65.35) and the marker
(23,900.57). It keys on σ alone, so it catches **both** populations — the zeroed means and
the real-mean-with-undetermined-direction cases.

That `"?"` reaches the database by two independent routes:

* canonical `reversibility` — `Estimate_Reaction_Reversibility.py:225` (this branch),
  during `./Estimate_Reaction_Reversibility.py EQ`
* the per-source operator in `thermodynamics['eQuilibrator'][2]` —
  `_thermo_helpers.py:414`, at energy-write time

**Result: 4,933 of 4,934 return `?`.** The single exception is `rxn11221`, an ABC
transporter, caught by `abc_transporter_heuristic` at cascade position 2. Its direction
comes from the sign of the ATP coefficient — a mechanistic call that never touches the
energy — so the bad number is irrelevant to it. Structural rules deliberately precede the
gate for exactly this reason; vetoing them on σ would discard a sound answer over a number
we are not using.

### `?` means eQuilibrator has no opinion, not that nobody does

All 4,933 retain their Group-Contribution and dGPredictor operators in `thermodynamics`,
untouched. Restoring coverage for any of them is a promotion-policy decision, not a
re-estimation — the material is already in the record.

---

## 5. Still open

1. **`is_using_group_contribution()`** — eQuilibrator exposes a boolean for exactly this
   condition. It would replace the σ threshold with eQuilibrator's own verdict and
   distinguish the two populations directly, which the σ test cannot.
2. **The retrieval script's MetaNetX collapse** — it keys the reaction formula on compound
   id and discards compartment, so 1,102 transport reactions (and 76 stereoisomer cases
   such as `rxn00816`, where D-glucose and galactose merge) carry a ΔG′° for a *different*
   reaction. Fixable at source.
3. **`standard_dg_prime_multi()`** — returns the full covariance that Beber 2022 introduced.
   The current per-reaction loop cannot; we store 25,028 independent scalar σ's, which the
   paper says systematically *overestimates* uncertainty.
4. **Canonical `deltag`** — still mostly the Group-Contribution value, so `reversibility`
   and `deltag` can differ in provenance. Reconciling them is
   `Promote_Reaction_Thermodynamics_to_Canonical.py`, deliberately not part of
   `Rerun_Thermodynamics.sh`.
