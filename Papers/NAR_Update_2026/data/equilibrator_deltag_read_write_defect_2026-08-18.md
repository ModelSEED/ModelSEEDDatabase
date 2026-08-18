# The `EQ` step stopped reading eQuilibrator energies in 2023

Fixed for `EQ` in [PR #285](https://github.com/ModelSEED/ModelSEEDDatabase/pull/285).
The same bug is still live in the `GC` path — see §5.
Code references are to `origin/dev` @ `49563c6f`.

**In one line:** a two-step handoff was broken in 2023 when someone changed the first
step and not the second, so the reversibility calculation has been scoring the wrong
method's energy ever since.

---

## 1. Every reaction stores energies in two places

```
rxn01442  ITP:D-glucosamine 6-phosphotransferase

  deltag        : -2.91          <- ONE canonical energy: "the official answer"
  reversibility : '='

  thermodynamics:                <- what each method SEPARATELY says
      Group contribution     : [-14.00, 8.85, '=']
      dGPredictor            : [ -2.91, 0.02, '=']
      dGPredictor-ModelSEED  : [ -2.73, 0.87, '=']
      eQuilibrator           : [ -2.31, 0.78, '=']
```

`thermodynamics` holds four independent opinions, side by side. `deltag` holds the single
canonical value the rest of the database uses.

Note that `deltag` here is `-2.91` — exactly **dGPredictor's** number. Not eQuilibrator's
`-2.31`. Hold onto that.

## 2. How it used to work (2019 → 2023)

Two steps, in order:

**Step A** — `Update_Reaction_eQuilibrator_Energies.py` wrote eQuilibrator's answer to
**both** places:

```
thermodynamics['eQuilibrator'] = -2.31     # its own slot
deltag                         = -2.31     # AND copied to canonical
```

**Step B** — `Estimate_Reaction_Reversibility.py EQ` read `deltag` to decide direction.

Step B read `deltag`, and Step A had just put eQuilibrator's number there. **So Step B was
reading eQuilibrator's energy.** Correct by construction.

## 3. What changed on 2023-09-13

Commit `3e50646b` made Step A write only its own slot:

```
thermodynamics['eQuilibrator'] = -2.31     # its own slot
deltag                         = (untouched)
```

**Step B was never updated.** It still reads `deltag`.

So `deltag` now holds whatever something else last put there — for `rxn01442`, dGPredictor's
`-2.91`. Step B reads that, derives a direction from it, and records it as the
eQuilibrator-informed answer.

> **Analogy.** Each method writes its result on its own sticky note, and there is a
> whiteboard holding "the official number". eQuilibrator used to put up its sticky note
> *and* copy the number to the whiteboard. The direction step reads the whiteboard. In
> 2023 eQuilibrator stopped copying to the whiteboard — but the direction step still reads
> it. So it reads whatever someone else left up there, and files the answer under
> eQuilibrator's name.

### Why nobody noticed for three years

1. **Nothing errors.** `deltag` always holds a valid number — just the wrong method's.
2. **The answers look reasonable.** `-2.91` vs `-2.31`: both plausible, usually the same
   direction.
3. **The label is always right.** The pipeline reports "eQuilibrator direction" regardless
   of whose number it used, so nothing downstream can tell.

### A correction to an earlier version of this note

It previously said the `EQ` step "never" read eQuilibrator energies. That was wrong: it
read them correctly for four years. The original design was sound; a later commit broke
it. The symptom today is the same, but the fix is a **restoration**, not a redesign.

---

## 4. The evidence

**Nothing in the pipeline writes `deltag` any more.** Scanning every Thermodynamics script,
the only writer is `Promote_Reaction_Thermodynamics_to_Canonical.py` — which is *not* in
`Rerun_Thermodynamics.sh` and is referenced by no script, wrapper, or CI config. It runs
only when typed by hand. `deltag` is therefore frozen data, not pipeline output.

**So `deltag` has drifted into a mix of sources and vintages.** Of 39,201 reactions with a
usable value:

| `deltag` matches | n | % |
|---|---:|---:|
| **no current source — stale** | **19,375** | **49.4%** |
| dGPredictor | 8,662 | 22.1% |
| Group contribution | 5,471 | 14.0% |
| dGPredictor-ModelSEED | 4,084 | 10.4% |
| eQuilibrator | 1,609 | 4.1% |

For the reactions the `EQ` step actually processes, `deltag` is eQuilibrator's number only
**7.2%** of the time (1,797 of 24,937).

**It changes answers.** Replaying both energies through the cascade, restricted to
reactions where eQuilibrator is confident (σ ≤ 5 kcal/mol), **1,575 directions differ**:

```
rxn00003  pyruvate:pyruvate acetaldehydetransferase (decarboxylating)
  canonical deltag       : 8.27 +/- 0.90     <- what the EQ step reads
  thermodynamics[eQuil.] : 7.14 +/- 0.89     <- what it should read
  dev answer  : '<'      true answer : '='
```

A 1.13 kcal/mol input difference flips reversible to irreversible-backward.

**Both updaters were changed in the same push**, two days apart, and neither reader was
updated — which is why the GC path has the identical bug:

| script | last wrote `deltag` | stopped by |
|---|---|---|
| `Update_Reaction_GroupContribution_Energies.py` | `17c9739b` (2023-09-11) | `6638af00`, 2023-09-13 |
| `Update_Reaction_eQuilibrator_Energies.py` | `17c9739b` (2023-09-11) | `3e50646b`, 2023-09-13 |

No later commit repaired either. `c887ede8` (2026-06-04) refactored around the defect;
`ad34d6ab` (2026-08-07) left it alone.

---

## 5. The same bug in the `GC` path — not fixed

`Estimate_Reaction_Reversibility.py GC` reads `deltag` too, and gets the right number only
**20.7%** of the time. But it fails *differently*, and the difference decides the fix.

Simulating the full pipeline both ways, **1,156 final directions change** — but only **5**
are hard flips (`>` ↔ `<`). The stale values here are mostly older GC numbers, a median of
2.07 kcal/mol from current ones, so the cascade usually lands the same way.

**The real cost is lost coverage: 843 reactions move from `?` to a confident call.**
`deltag` is the `10000000` sentinel for **892 reactions that have a perfectly good stored
GC energy**, so the step reports "no data" and returns `?`:

```
rxn31387  R_NADH_monodehydroascorbate_oxidoreductase
   deltag      = 10000000.0 (sentinel)  -> current answer '?'
   GC sublist  = [-2.57, 12.67]         -> correct answer '='
```

### Why the GC fix is not the same one-line swap

Of those 892, three carry `|dG| > 1000` kcal/mol and five carry error > 100 —
`rxn31466` sits at **−3505.99 ± 244.6**. They are invisible today *because*
`Promote_Reaction_Thermodynamics_to_Canonical.py` applies `MAX_ABS_DG = 1000` and
`MAX_ERR = 100` before writing `deltag`.

> **`deltag` is accidentally acting as a quality filter.** Reading it is wrong, but
> protective. Swapping the GC path to the per-source value without porting those guards
> would let all 892 through, nonsense included.

The EQ fix could be a clean swap because the EQ rule set carries its own quality gate
(`eq_undecomposable_heuristic`). The GC cascade has no equivalent — `stored_bounds_heuristic`
abstains on a wide error bar rather than rejecting it.

A follow-up PR should: swap the GC source, add a plausibility gate to the GC rule set
mirroring `Promote_*`, and report the ~884 recovered reactions as an intended coverage
gain. Kept out of PR #285 so the byte-identical GC column stays a control proving the
7,103 eQuilibrator changes come from the new rule set and not the refactor.

`DGP` shares the pattern but is inert — that step is not in `Rerun_Thermodynamics.sh`.

---

## 6. The decision this needs

`deltag` is currently **neither** pipeline output nor a curated field — a frozen artifact
that live code reads as authoritative. Two coherent resolutions:

1. **Make it output again** — add `Promote_Reaction_Thermodynamics_to_Canonical.py` to
   `Rerun_Thermodynamics.sh`, so `deltag` is regenerated from a stated policy every run.
2. **Declare it curated** — document that it is set only by explicit promotion, and convert
   the remaining readers to per-source sources so no estimation step depends on it.

Either closes the gap. Doing neither means the drift keeps widening with every GC rebuild.
This is a policy question about what `deltag` is *for*, so it is not one I have decided.

---

## Summary

| | |
|---|---|
| **What** | The `EQ` reversibility step scores canonical `deltag`, not eQuilibrator's own estimate |
| **Since** | `3e50646b`, 2023-09-13 — the writer stopped setting `deltag`; the reader was never updated |
| **Scale** | `deltag` is eQuilibrator's number 7.2% of the time; ≥1,575 directions differ |
| **Cause** | A two-step handoff where only the first step was changed |
| **Fixed** | `EQ` path, PR #285 |
| **Open** | `GC` path — same bug, 892 reactions with an unreachable energy, needs a plausibility gate first (§5); and `deltag`'s role needs a decision (§6) |

Companions: `equilibrator_reversibility_heuristics_2026-08-18.md` (the fix as built),
`equilibrator_zero_energy_defect_2026-08-18.md` (the separate zeroed-energy defect).
