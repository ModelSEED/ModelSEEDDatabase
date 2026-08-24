# The reversibility index, and what it tells us about ModelSEED's directions

**2026-08-21.** Data: ModelSEED dev @ `49563c6f`. Branch `eQuilibrator-fix`.
Companions: `equilibrator_reversibility_heuristics_2026-08-18.md`,
`equilibrator_zero_energy_defect_2026-08-18.md`.

---

## The short version

We rebuilt eQuilibrator's 2012 reversibility index as a standalone function and
used it to ask two questions: *are our structural shortcuts doing real work?* and
*how much of a direction call is the rule set rather than the chemistry?*

- **The index cannot replace the ATP-synthase and ABC-transporter shortcuts.**
  It looks like it agrees with them; the agreement is an artifact.
- **On eQuilibrator energies the old and new rules agree only 71.6%** — almost
  all of the gap is reactions eQuilibrator never actually scored.
- **The four sources look more consistent than they are.** Under the current
  defaults most of their agreement is both sides saying "reversible".
- **On dGPredictor they agree 85–90%, for a bad reason:** nearly half of the old
  method's directional calls on the high-uncertainty source come from rules that
  never look at the error bar.
- **eQuilibrator's way of saying "I don't know" is invisible** to anything
  looking for ModelSEED's. We now translate between them.

Nothing was applied to the live database — every number here is a dry run.

---

## 1. What the reversibility index is

Γ answers one question: **by what factor would every reactant concentration have
to change before the reaction ran backwards?** More than 1000× and physiology
can't realistically do it, so we call the reaction irreversible.

```
ln Γ = (2 / Σ|ν|) · ΔG′m / RT
```

ΔG′m is the free energy with every reactant at 1 mM. Water and protons are left
out — the cell buffers both. The `2/Σ|ν|` factor makes a two-substrate reaction
comparable to a five-substrate one.

It now lives in `Scripts/Thermodynamics/reversibility_index.py` as a plain
function. Checked against eQuilibrator's own published numbers: **17,776
reactions agree.**

**One rule throughout: a transport reaction is just a reaction.** A metabolite in
two compartments is two species — no membrane terms, no special cases. That is
what lets transport and cytosolic reactions be compared on the same axis.

---

## 2. The index does not rescue ATP synthase or ABC transporters

ModelSEED decides these two families *structurally*, by pattern-matching the
reaction before looking at any energy. That has always felt like a hack, so we
scored them with the index instead to see whether thermodynamics could do the job.

**ATP synthase (15 reactions).** Every one comes out at `|ln Γ| = 11.97` — not
similar, *identical*. The translocated protons cancelled when the formula was
built, so eQuilibrator scored `ADP + Pi ⇌ ATP + H₂O` fifteen times. The index
calls the cell's most famously reversible enzyme irreversible, splitting 7
forward / 8 backward purely on **how each reaction was written down**. (Letting
the protons count returns 13 of 15 to "reversible".)

**ABC transporters (1,150 reactions).** 94% agreement with the structural rule,
which sounds like success. It isn't: the transported substrate also cancels,
leaving plain ATP hydrolysis at `ln Γ = −7.18` against a cut-off of −6.91. It
clears the bar by 4%, and **1,074 of the 1,111 hard calls sit within 15% of the
threshold** — one rounding from falling apart.

**Conclusion:** the shortcuts exist because the energy behind these reactions
can't be trusted. The index doesn't fix that; it makes it measurable.

---

## 3. Rule set vs. chemistry

`Compare_Reversibility_Heuristics.py` hands two rule sets the *exact same*
energy, so any disagreement is the rules and nothing else. (Comparing the shipped
GC directions against the shipped eQuilibrator directions can't tell you this —
those runs differ in their inputs too.)

### The two methods ask different questions

**GC asks a yes/no question about a range.** It computes ΔG at the extremes of a
concentration box (10 µM–20 mM per reagent), widens it by the error bar, and asks
*does the sign hold everywhere?* If not it falls through to a band test, then to
a rule scoring **chemistry rather than energy** (does the reaction consume ATP,
release CO₂?), then to a bare default of "reversible".

**The index asks a how-much question about a point.** Everything at 1 mM, one
continuous number, past 1000× it's irreversible. No fallbacks, no chemistry rule.

So the GC cascade can reach an answer without ever consulting the
thermodynamics, and the index cannot. That turns out to matter.

### On eQuilibrator energies

The two rule sets **agree on 71.6%** of 25,028 reactions, and the gap is one cell:

> **4,369 reactions the old rules call reversible (`=`), the new rules call
> unknown (`?`).**

Those are reactions eQuilibrator declined to estimate. The old rules never looked
at the error bar, so they sailed past.

### The confidence margin

The plain 2012 rule is `|ln Γ| > ln(1000)`. The margin adds one term:

```
|ln Γ| − z·σ(ln Γ) > ln(1000)
```

In words: **the best guess clearing the 1000× bar isn't enough — the pessimistic
end of the error bar has to clear it too.** `z = 0` is Noor 2012 as published;
`z = 1` demands one sigma of headroom, in the spirit of eQuilibrator 3.0.
`σ(ln Γ)` is just the reported energy uncertainty pushed through the same
`2/Σ|ν|/RT` scaling as the energy itself.

It only bites when that uncertainty is large next to the threshold it must clear:

| source | median σ(ln Γ) | vs the ln(1000) threshold | reactions whose 1σ exceeds the whole threshold |
|---|---:|---:|---:|
| `dGPredictor` | 0.30 | 0.04× | 0% |
| `eQuilibrator` (decomposable) | 0.57 | 0.08× | 0% |
| `dGPredictor-ModelSEED` | 18.43 | **2.67×** | **70%** |

That last row is the whole story. A worked example — `rxn00985`, propionate
kinase, on dGPredictor-ModelSEED:

```
ΔG′°     =  9.25 ± 11.82 kcal/mol      Σ|ν| = 4
ln Γ     =  (2/4) ×  9.25 / RT  =  7.81
σ(ln Γ)  =  (2/4) × 11.82 / RT  =  9.98

point estimate:   7.81         > 6.91   ->  "<"   (irreversible)
with one sigma:   7.81 − 9.98  = −2.17  ->  "="   (no call)
```

The best guess says this runs backwards; the error bar says ΔG could just as
easily be negative. Same record, two defensible readings.

**"=" therefore covers two situations**, and the rule set records which in the
status string so they never get conflated: `RI:reversible` (confidently inside
the window) vs `RI:ambiguous` (the interval straddles the threshold). On
dGPredictor-ModelSEED that split is **5,687 confident vs 19,512 ambiguous**
(plus 76 with nothing left to count). A model reading the database sees a bare
`=` for all 25,275 and cannot tell them apart.

### On dGPredictor energies

dGPredictor had been falling through to the GC rules by default, so this is the
comparison that changes what ships. Against the index with the margin (`DGP`,
z = 1) and without it (`RI`, z = 0):

| | reactions | GC vs `DGP` | GC vs `RI` |
|---|---:|---:|---:|
| `dGPredictor` (σ median 0.35) | 27,715 | **89.8%** agree | 89.1% agree |
| `dGPredictor-ModelSEED` (σ median 21.17) | 31,924 | **84.7%** agree | 62.7% agree |

**On plain dGPredictor the margin is irrelevant** — error bars are small enough
that it changes almost nothing. Of the 2,816 reactions that differ, the dominant
cell is 2,078 that GC calls forward and the index calls reversible.

**On dGPredictor-ModelSEED the margin does all the work, and the comparison
inverts.** With it, the gap is 4,457 reactions GC calls forward and the index
won't. Without it, the gap flips to 9,291 that GC calls reversible and the index
calls directional — now the index is the bolder one. How many hard directional
calls each produces:

| | GC | index + margin | index, point estimate |
|---|---:|---:|---:|
| `dGPredictor` | 58.7% | 52.4% | 53.6% |
| `dGPredictor-ModelSEED` | 34.9% | 20.8% | 56.0% |

**GC sits in the middle — but not because it is better calibrated.** It mixes
rules that respect the error bar with rules that ignore it. Of the 9,922
reactions it calls forward on dGPredictor-ModelSEED:

- 5,389 (54%) from the concentration-bounds test, which *does* use σ;
- 3,267 (33%) from the low-energy-compound rule, which **never looks at σ**;
- 1,266 (13%) from the ABC shortcut, which never looks at the **energy** at all.

The shift between sources shows the mechanism. Moving from the low-σ source to
the high-σ one, GC's bounds test correctly backs off — from deciding 49% of
reactions to 20%. But the reactions it declines fall through to the σ-blind rule,
which more than doubles its forward calls (1,311 → 3,267) as the evidence gets
*weaker*. GC responds to the error bars in one rule and partly undoes it in the
next.

**Transport is the other asymmetry**, because GC's ABC shortcut fires before any
energy is read and the index has no equivalent:

| | non-transport | transport |
|---|---:|---:|
| `dGPredictor` | 91.9% agree | **77.4%** agree |
| `dGPredictor-ModelSEED` | 86.8% agree | **73.8%** agree |

---

## 4. Do the four sources agree with each other?

The comparisons so far hold the energy fixed and vary the rules. This one does
the opposite: **hold the rules fixed and vary which prediction supplies the
energy.** Each source is scored from its own `thermodynamics[source]` entry, and
each cell is the percentage of shared reactions where two sources give the same
direction.

Two configurations. First, what ModelSEED does today — the default cascade on
all four. Then the realistic alternative: **Group contribution keeps the cascade
it was designed around** (it is also the byte-compare anchor), and only the three
prediction sources move to the index.

`./Compare_Reversibility_Heuristics.py --matrix --sets GC RI`

**A. Default cascade on all four — % identical direction**

| | GroupContrib | eQuilibrator | dGPredictor | dGP-ModelSEED |
|---|---:|---:|---:|---:|
| **GroupContrib** | — | 75.7% | 62.2% | **83.3%** |
| **eQuilibrator** | 75.7% | — | 66.7% | 76.0% |
| **dGPredictor** | 62.2% | 66.7% | — | 62.2% |
| **dGP-ModelSEED** | **83.3%** | 76.0% | 62.2% | — |

**B. Index on the three prediction sources, GC cascade on Group contribution**

| | GroupContrib | eQuilibrator | dGPredictor | dGP-ModelSEED |
|---|---:|---:|---:|---:|
| **GroupContrib** | *(GC rules)* | 62.2% | 54.3% | 61.6% |
| **eQuilibrator** | 62.2% | — | 51.6% | 55.0% |
| **dGPredictor** | 54.3% | 51.6% | — | 59.2% |
| **dGP-ModelSEED** | 61.6% | 55.0% | 59.2% | — |

Overlaps range from 17,021 (eQuilibrator ∩ dGPredictor) to 25,969
(GroupContrib ∩ dGP-ModelSEED); no cell rests on fewer than 17,000 reactions.

### Agreement falls everywhere, and it is abstention that disappears

Average pairwise agreement drops from **71.0% to 57.3%**, and every single pair
moves the same way. Splitting each pair's agreement into *why* they agree shows
what is actually being lost:

| pair | agreement | both said `=` | direct conflict |
|---|---:|---:|---:|
| GroupContrib / eQuilibrator | 75.7 → 62.2% | 39.6 → 32.4% | 0.6 → 0.7% |
| GroupContrib / dGPredictor | 62.2 → 54.3% | 33.7 → 32.6% | 2.7 → 3.2% |
| GroupContrib / dGP-ModelSEED | **83.3 → 61.6%** | **54.1 → 34.2%** | 0.4 → 2.2% |
| eQuilibrator / dGPredictor | 66.7 → 51.6% | 31.7 → 27.4% | 3.3 → 4.0% |
| eQuilibrator / dGP-ModelSEED | 76.0 → 55.0% | 43.9 → 24.8% | 0.2 → 3.0% |
| dGPredictor / dGP-ModelSEED | 62.2 → 59.2% | 34.5 → 27.9% | 1.9 → **5.9%** |

The best-agreeing pair in the default matrix is Group contribution and
dGP-ModelSEED at 83.3% — and **54.1 of those 83.3 points are both sources saying
"reversible".** They agree by both declining to answer; only 29.2% is two sources
committing to the same direction.

That is the pattern throughout. In every pair the mutual-`=` share falls and
direct conflicts — one source says forward, the other says backward — rise, in
one case from 0.2% to 3.0%. **The conflicts were always there; the default
cascade's fallback to "reversible" was hiding them.**

Worth being precise about what does *not* happen: agreed-upon directions don't
rise to compensate. They fall in five of the six pairs (only dGPredictor /
dGP-ModelSEED gains, 27.7 → 31.3%). The index isn't converting abstention into
consensus — it is converting it into a mix of consensus and conflict, because
the two rule families commit on *different* reactions.

**So the current defaults make the four thermodynamic sources look more
consistent than they are.** Much of that consistency is silence, and the pair
that looks most alike is the one that abstains most.

### One caveat on the eQuilibrator row

Under the index, eQuilibrator returns `?` for the 4,934 reactions it declined
(see §5), and a `?` matches nothing. That drags its row down for a reason that
isn't disagreement. Restricted to reactions where both sources actually made a
call:

| pair | as tabulated | excluding `?` |
|---|---:|---:|
| GroupContrib / eQuilibrator | 62.2% | **72.9%** |
| eQuilibrator / dGPredictor | 51.6% | **63.3%** |
| eQuilibrator / dGP-ModelSEED | 55.0% | **68.4%** |

The three pairs involving no `?` still fall, so the finding holds — but
eQuilibrator is a better-behaved source than its raw row suggests, once an honest
refusal stops counting as a wrong answer.

---

## 5. Two ways of saying "I don't know"

| | how it reports "no value" |
|---|---|
| **ModelSEED** | `[10000000.0, 10000000.0, '?']` — impossible to mistake |
| **eQuilibrator** | a normal-looking signed energy, with σ = 1e5 kJ/mol |

eQuilibrator never raises an error. When it can't decompose a reaction it returns
**zero**, discards the computed value, then — in the same call — adds the pH
correction on top anyway. What reaches us is `0 + pH correction`: a plausible,
signed, nonzero number whose entire content is how many protons the reaction
releases at pH 7.

The only reliable signature is the error bar. The largest *genuine* one in the
table is 65 kcal/mol against a marker of 23,900, and **nothing at all sits in
between** — a three-order-of-magnitude gap with no borderline case to get wrong.

Two scripts came out of this:

- **`Check_eQuilibrator_Energy_Errors.py`** finds them: **4,934 declined
  reactions**, 632 with a bare zero energy, and 1,178 where our own retrieval
  step scored the wrong reaction. Self-tested.
- **`Normalize_eQuilibrator_Sentinels.py`** translates eQuilibrator's marker into
  ModelSEED's, so anything that already skips one will skip the other. Dry-run by
  default.

### The judgment call

288 records have an energy of exactly `0.0` with a *normal* error bar. Tempting
to blank those too — but they aren't the failure branch, which always stamps the
1e5 marker. 280 have an error bar of exactly zero and 232 are simple
isomerizations (L-lysine ⇌ D-lysine, 16α- ⇌ 16β-hydroxysteroid). **eQuilibrator
can't see stereochemistry**, so both sides decompose to the same groups and the
difference is exactly zero with exactly zero error.

That is a real statement — these pairs *are* isoenergetic — not a failure report.
So we leave them alone by default.

---

## What changed in existing code

- `reversibility_heuristics.py` — added the index-based rule sets; `Context` now
  delegates to the shared function so the two can't drift apart.
- Both dGPredictor sources use the index instead of the GC fallback.
- `test_reaction_direction.py` — **Group contribution is now the only
  byte-compare anchor.** dGPredictor was one too, until we changed it on purpose.
  It still passes exactly.
- Everything reads energies from the per-source `thermodynamics` dictionary,
  never the flat `deltag` field — nothing has rewritten that field since the
  additive refactor, so it holds whichever source last promoted to canonical.

## What we didn't do

- **Nothing was applied to the live database.** `--apply` was verified end-to-end
  on a sandbox copy, where it changed only the eQuilibrator entries and left
  everything else byte-identical.
- The 1,178 wrong-formula reactions are reported, not repaired. That needs the
  reaction formula rebuilt with compartments intact.
- 4 of the 4,934 declined reactions are also the current canonical `deltag`.
  Flagged, not touched; picking a replacement belongs to
  `Promote_Reaction_Thermodynamics_to_Canonical.py`.
