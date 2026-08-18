# The `EQ` reversibility step is not reading eQuilibrator energies

**A broken read/write contract in the Thermodynamics pipeline.**
Introduced 2023-09-13, live on `dev` continuously since. Fixed for `EQ` in
[PR #285](https://github.com/ModelSEED/ModelSEEDDatabase/pull/285); **the same pattern
remains in the `GC` and `DGP` paths and is analysed in §3.**

Code references are to `origin/dev` @ `49563c6f`.
Companions: `equilibrator_reversibility_heuristics_2026-08-18.md` (the fix),
`equilibrator_zero_energy_defect_2026-08-18.md` (the separate zeroed-energy defect).

---

## 1. Verifying the problem

Four checks, each independently reproducible against `origin/dev`.

### Check 1 — Nothing in the pipeline writes canonical `deltag`

Scan every Thermodynamics script for a canonical write:

```
$ for f in Scripts/Thermodynamics/*.py; do
    grep -cE "\['deltag'\]\s*=|\['deltagerr'\]\s*=" $f
  done
```

```
Promote_Reaction_Thermodynamics_to_Canonical.py : 2 writes
(no other script writes it)
```

`Rerun_Thermodynamics.sh` runs nine scripts. **`Promote_…` is not one of them.**

> **`deltag` is not pipeline output.** Running the full pipeline today cannot change it.
> It is static data, frozen from whenever `Promote_…` was last invoked by hand.

### Check 2 — But the `EQ` step reads `deltag` as its energy source

`Estimate_Reaction_Reversibility.py:74-75`:

```python
    status, thermoreversibility, source_label = run_reversibility(
        rxn_entry, top_level_energy(db_level), DEFAULT_HEURISTICS)
```

resolving to `reversibility_heuristics.py:75-95`:

```python
def _energy_for(rxn_entry, db_level):
    rxn_dg = rxn_entry['deltag']                      # L80  <-- the canonical field
    ...
    if db_level:
        if not _is_source_eligible(rxn_entry, db_level):
            return None, None, None                   # db_level only GATES
        label = DB_LEVEL_LABEL[db_level]
        append = label if _thermo_pair(rxn_entry, label) is not None else None
        return rxn_dg, rxn_dge, append                # L95  <-- deltag, tagged as `label`
```

Line 94 looks up `_thermo_pair(rxn_entry, label)` — **it has the eQuilibrator value in
hand** — then discards it and returns `rxn_dg` from `deltag`, labelled eQuilibrator.
`db_level` selects *which reactions are processed*, never *which energy is used*.

### Check 3 — So what is actually in `deltag`?

Matching every reaction's canonical `deltag` against each stored per-source value:

| `deltag` exactly matches | n | % |
|---|---:|---:|
| dGPredictor | 8,662 | 22.1% |
| Group contribution | 5,471 | 14.0% |
| dGPredictor-ModelSEED | 4,084 | 10.4% |
| eQuilibrator | 1,609 | 4.1% |
| **none of them** | **19,375** | **49.4%** |

*(39,201 reactions with a usable `deltag`.)*

Half of all canonical energies correspond to **no value currently in the database**. For
those, the distance to the *current* Group-Contribution value is a median of
2.07 kcal/mol — consistent with them being pre-Convention-A GC numbers, stranded when
`ad34d6ab` rebuilt GC on 2026-08-07 without re-promoting.

Restricting to the reactions the `EQ` step actually processes:

| EQ-eligible reactions with a usable `deltag` | 24,937 |
|---|---:|
| `deltag` == eQuilibrator value **(correct)** | **1,797 — 7.2%** |
| `deltag` == Group-contribution value | 1,446 — 5.8% |
| `deltag` == neither | 21,694 — 87.0% |

**92.8% of the time, the `EQ` run scores a number that is not eQuilibrator's.**

Corroborating: canonical `deltagerr` is **never** the eQuilibrator undecomposable
sentinel — 0 cases out of 4,934 reactions whose eQuilibrator record carries it. The
canonical error bars did not come from eQuilibrator either.

### Check 4 — It changes the answer

Replaying dev's `EQ` logic against `deltag`, versus the same cascade against
eQuilibrator's own energy, on reactions where eQuilibrator is confident (σ ≤ 5 kcal/mol):
**1,575 reactions get a different direction.**

```
rxn00003  pyruvate:pyruvate acetaldehydetransferase (decarboxylating)
          (1) CO2 + (1) ALCTT <= (2) Pyruvate + (1) H+

  canonical deltag       :   8.27 +/- 0.90     <- what dev's EQ run reads
  thermodynamics[eQuil.] :   7.14 +/- 0.89     <- what it should read

  dev  EQ answer : '<'   via  lowE: 9.63:1.0
  true EQ answer : '='   via  EQ:reversible: 6.03 +/- 0.75
```

A 1.13 kcal/mol difference in input flips a reversible reaction to irreversible-backward.
`rxn00010` and `rxn00011` behave the same way. The errors are small in energy and
decisive in outcome, because the cascade thresholds sit close to these values.

**Conclusion: verified.** The `EQ` reversibility step reads a stale, heterogeneous
canonical field instead of the eQuilibrator estimate, and it changes at least 1,575
directions among confidently-estimated reactions alone.

---

## 2. Where the bug came from

### The original design was correct

Until 2023-09-13 a two-sided contract held.

**Writer** — `Update_Reaction_eQuilibrator_Energies.py` @ `17c9739b` (2023-09-11):

```python
    if(float(eq_reactions[rxn]['dge']) > 100):
        continue                                                       # quality guard

    reactions_dict[rxn]['deltag']=float(eq_reactions[rxn]['dg'])       # establishes it
    reactions_dict[rxn]['deltagerr']=float(eq_reactions[rxn]['dge'])
    if('EQU' not in notes_list):
        notes_list.append('EQU')                                       # marks it
```

**Reader** — `reversibility_heuristics.py` (unchanged, then and now): read `deltag`, gate
on `EQU`.

The contract: *`deltag` holds the eQuilibrator energy wherever `EQU` is set.* The reader
was correct **because** the writer guaranteed it. Reading `deltag` was reading
eQuilibrator energies, exactly as intended.

### `3e50646b` (2023-09-13) deleted the writer's half

```python
    if(float(eq_reactions[rxn][1]) > 100):
        pass                                                # guard -> no-op

    reactions_dict[rxn]['thermodynamics'][label]=eq_reactions[rxn]
                                                            # deltag write DELETED
                                                            # EQU stamp DELETED
```

The reader was not touched. Its precondition is now established by **nobody**.

### Why it failed silently

1. **`EQU` is persisted data, not recomputed state.** The 17,094 notes written before
   2023-09-13 remain in the JSON forever, so the gate keeps passing. Only **1,586** of
   them still have `deltag` equal to the eQuilibrator value.
2. **`deltag` was never emptied.** It still holds a valid float, so `_energy_for` returns
   a number instead of `(None, None, None)`.
3. **The cascade cannot detect provenance.** `stored_bounds`, `mMdeltaG`, `lowE` are
   arithmetic on a float. Nothing checks where it came from.

Note `EQU` ≠ `EQC`. `EQU` was written by the *updater* to mean "this deltag is
eQuilibrator's". `EQC` is written by the *retrieval* script to mean "all reagents had
eQuilibrator structures" — a different claim, still maintained (27,024 reactions).

### It was never repaired

| commit | date | writes `deltag`? | stamps `EQU`? | guard |
|---|---|---|---|---|
| `17c9739b` | 2023-09-11 | yes | yes | `continue` |
| **`3e50646b`** | **2023-09-13** | **no** | **no** | `pass` (dead) |
| `c887ede8` | 2026-06-04 | no | no | none |
| `ad34d6ab` | 2026-08-07 | no | no | none |

`c887ede8` ("Add additive per-source thermodynamic descriptions") refactored *around* the
defect; `ad34d6ab` left it alone. Live continuously from **2023-09-13 to 2026-08-18** —
just under three years — reaching `dev` via the merges of 2024-01-16, 2024-05-02 and
2024-12-11.

**The drift grows over time.** Every unrelated event that touches canonical energies —
the PR #265 promotion, the Convention A GC rebuild — moves `deltag` further from the
eQuilibrator value while the `EQ` step keeps reading it. That is why only 7.2% still
agree.

---

## 3. Where it persists in the current code

### Fixed on `eQuilibrator-fix`: the `EQ` path

`reversibility_heuristics.py` now routes the EQ run to the per-source value:

```python
def energy_source_for_level(db_level):
    if db_level == 'EQ':
        return per_source_energy(DB_LEVEL_LABEL['EQ'])   # thermodynamics['eQuilibrator']
    return top_level_energy(db_level)                    # everything else: unchanged
```

Worth noting: **`per_source_energy` already existed on `dev`** at
`reversibility_heuristics.py:285`, fully written and never called — its only appearances
are the module docstring and a re-export. The building block was there; nothing was wired
to it.

Side effect, and a desirable one: `EQU` stops being consulted. Eligibility becomes "does
`thermodynamics['eQuilibrator']` hold a non-sentinel pair" — live data rather than a 2023
fossil.

### NOT fixed: the `GC` path has the identical defect

`top_level_energy` is still the source for every other level, so
`Estimate_Reaction_Reversibility.py GC` reads `deltag` too — same line, same mechanism:

| GC-eligible reactions with a usable `deltag` | 26,421 |
|---|---:|
| `deltag` == Group-contribution value **(correct)** | **5,471 — 20.7%** |
| `deltag` == eQuilibrator value | 1,481 — 5.6% |
| `deltag` == neither | 19,469 — 73.7% |

**The GC run reads the wrong energy 79.3% of the time.**

#### But it fails differently, and that matters

Simulating the full pipeline both ways — GC step, then EQ step, in
`Rerun_Thermodynamics.sh` order — and comparing final canonical `reversibility`:

| transition | n |
|---|---:|
| `?` → `=` | 608 |
| `?` → `>` | 213 |
| `>` → `=` | 105 |
| `<` → `=` | 43 |
| `=` → `?` | 43 |
| `=` → `>` | 41 |
| `>` → `?` | 32 |
| `<` → `?` | 30 |
| `?` → `<` | 22 |
| `=` → `<` | 14 |
| `>` → `<` | 4 |
| `<` → `>` | 1 |
| **total** | **1,156** |

Only **5** are hard directional flips (`>` ↔ `<`), and 105 turn a confident call into
`?`. So the eQuilibrator failure mode — confident wrong answers — is largely *absent*
here, because the stale values are mostly older GC numbers sitting a median of
2.07 kcal/mol from the current ones. Close enough that the cascade usually lands the same
way.

**The dominant effect is lost coverage: 843 reactions move from `?` to a confident
call.** The cause is a different branch of the same defect — `deltag` is the sentinel
`10000000` for reactions that nonetheless have a perfectly good stored GC energy, so
`_energy_for` returns `(None, None, None)`, `_incomplete_decision` fires, and the answer
is `?`:

```
rxn31387  R_NADH_monodehydroascorbate_oxidoreductase
   canonical deltag = 10000000.0  (sentinel)   -> current answer '?'
   GC sublist       = [-2.57, 12.67]           -> correct answer '='
```

**892 reactions** have a stored Group-Contribution energy that is unreachable through
`deltag`.

#### Why the GC fix is not a one-line swap

Those 892 are not all worth recovering:

| | n |
|---|---:|
| \|dG\| > 1000 kcal/mol — implausible | 3 |
| error > 100 kcal/mol — uselessly vague | 5 |
| would pass `Promote_*`'s guards | 884 |

`rxn31466` carries a GC energy of **−3505.99 ± 244.6 kcal/mol**. It is currently
invisible precisely *because* `Promote_Reaction_Thermodynamics_to_Canonical.py` applies
`MAX_ABS_DG = 1000.0` and `MAX_ERR = 100.0` before writing `deltag`.

> **`deltag` is accidentally functioning as a quality filter.** Reading it is wrong, but
> protective. Swapping the GC path to `per_source_energy` without porting those
> plausibility guards would let all 892 through, including the nonsense.

That is why the eQuilibrator fix could be a clean swap and the GC fix cannot: the EQ rule
set carries its own quality gate (`eq_undecomposable_heuristic`), so it does not depend on
`deltag` for filtering. The GC cascade has no equivalent — its `stored_bounds_heuristic`
simply abstains on a wide error bar rather than rejecting it.

#### Why it is out of scope for PR #285

Changing the GC path in the same PR would destroy the control that makes the eQuilibrator
change reviewable. As submitted, the GC and dGPredictor operator columns are
**byte-identical** to `dev`, which is what demonstrates the 7,103 eQuilibrator changes
come from the new rule set rather than from the refactor. Mixing in ~1,156 GC changes
would remove that evidence.

A follow-up PR should: swap the GC path to `per_source_energy('Group contribution')`, add
a plausibility gate to the GC rule set mirroring `Promote_*`'s `MAX_ABS_DG`/`MAX_ERR`, and
report the ~884 recovered reactions as an intended coverage gain.

#### `DGP`

Shares the pattern, but `Estimate_Reaction_Reversibility.py DGP` is not in
`Rerun_Thermodynamics.sh`, so it is currently inert. It should be fixed at the same time
as GC to stop the pattern reappearing if that step is ever enabled.

### The underlying ambiguity, for a decision

`deltag` is currently **neither** pipeline output nor a curated field — it is a frozen
artifact that live code still reads as though it were authoritative. Two coherent
resolutions:

1. **Make it output again** — add `Promote_Reaction_Thermodynamics_to_Canonical.py` to
   `Rerun_Thermodynamics.sh` so `deltag` is regenerated from a stated policy on every run.
2. **Declare it curated** — document that `deltag` is set only by explicit promotion, and
   convert the remaining `top_level_energy` readers to per-source sources so no estimation
   step depends on it.

Either closes the gap. Doing neither means the drift keeps widening. This is a policy call
about what `deltag` is *for*, so I have not made it.

---

## 4. Summary for review

| | |
|---|---|
| **What** | The `EQ` reversibility step scores canonical `deltag`, not the eQuilibrator estimate |
| **Since** | `3e50646b`, 2023-09-13 — the writer stopped setting `deltag`; the reader was never updated |
| **Scale** | 92.8% of EQ-eligible reactions score a non-eQuilibrator number; ≥1,575 directions differ |
| **Root cause** | A two-sided contract where only one side was changed, with no mechanism to detect the breach |
| **Fixed** | `EQ` path, in PR #285 |
| **Outstanding** | `GC` path, same defect: 79.3% wrong energy, 1,156 directions affected, 892 reactions with an unreachable GC energy. Needs a plausibility gate, not just a swap — see §3 |
| **Also outstanding** | `deltag`'s role is undefined: neither pipeline output nor a curated field |
