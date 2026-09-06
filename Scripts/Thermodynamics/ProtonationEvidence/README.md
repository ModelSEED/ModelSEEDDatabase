# Protonation evidence

How well founded is the acid-dissociation model behind each reaction's
thermodynamics, and which reactions should not be read at face value.

eQuilibrator's Legendre transform needs a **macroscopic** pKa ladder per
compound: an ordered list in which consecutive species differ by exactly one
proton. The count of entries above the reported pH sets the proton count of the
major microspecies, so the *composition* of the ladder changes the transformed
energy directly, not merely its precision. Getting the ladder wrong is not a
precision problem; it moves the answer.

These scripts establish where each compound's ladder comes from, whether that
provenance is sound, and what it implies for every reaction the compound
appears in.

## Scripts

| script | what it answers |
|---|---|
| `check_mobile_h.py` | Which of our structures have a ladder available in the upstream cache that InChIKey matching misses, and is taking it safe? |
| `build_ladder_table.py` | What macroscopic ladders can be cited from the literature, and what shape does the magnesium guard demand? |
| `validate_ladders.py` | Does every literature row carry a citation, and would each ladder actually be installed? |
| `grade_protonation_evidence.py` | Gold / silver / bronze per reaction. |

`paths.py` resolves everything. In-repo paths derive from the file location; the
eQuilibrator tree, its caches and the pKa-experiment tables come from
environment variables with documented defaults, because they live on the
analysis host. Missing inputs raise rather than being skipped — a grading run
that silently dropped half its inputs would report a flattering distribution.

## Why structural matching, and why it is fussy

The upstream cache is keyed by InChIKey, and exact-key matching misses ladders
that are genuinely present. FAD is the example: our InChI and the cache's agree
on formula and connectivity character-for-character, but we keep the
exchangeable protons in four separated groups where the cache merges all eleven
into one, which changes the first-block hash. Nothing about either structure is
wrong.

Matching modulo the `/h` layer would recover these — but `/h` also encodes
tautomers, and a tautomer can have genuinely different pKas. `check_mobile_h.py`
therefore compares the fixed-H spec, the mobile proton count, and the set of
atoms those protons may occupy, and accepts only a regrouping of the same
protons over the same atoms.

**Candidate selection must search, not take the first hit.** IPP and DMAPP are
constitutional isomers that InChI separates *only* in `/h` — same formula, same
connectivity, the double-bond position expressed as which atoms carry hydrogen —
and the cache holds a row for each. Taking the first candidate on a
formula+connectivity match compares DMAPP against the IPP row and rejects a
mismatch that selection itself manufactured. 79 compounds have more than one
candidate on the same formula, connectivity and stereo.

Result over the polyprotic-collapse population: 254 compounds safe with `/h`
identical, 12 safe as regroupings, **105 refused** on a fixed-H mismatch. The
refused set is chemically coherent — anthocyanins, curcuminoids and two CoA
thioesters, where flavylium/quinoidal and keto/enol assignment genuinely
differs.

## Grading

Reaction grade is the **minimum over its compounds**, which also captures the
mixed-enumeration defect: a reaction with one cache-enumerated and one
predictor-enumerated partner grades down, and that inconsistency is what moved
`rxn08789` by 355 kcal/mol. 94.3% of reactions mix grades.

Two schemes are emitted per reaction, because they answer different questions.

**Provenance** (`grade` column) — where the ladder came from. Honest, and a poor
discriminator: 1.3% gold, 79.5% bronze, because almost every reaction touches at
least one predicted compound.

**Mechanism** (`mechanism_grade` column) — can the uncertainty reach the pH 6–8
window, where a differing pKa moves dG'° roughly 8x more than one far from
neutrality? Gold 20.5%, silver 36.5%, bronze 43.0%.

A compound showing the symmetry-collapse signature counts as near-neutral
**whatever its predicted values**, because collapse means the real dissociations
are unknown: the predictor returns 2.11 for orthophosphate, whose true pKa2 is
7.20. Skipping that correction grades 2,843 compounds and 10,233 reactions too
generously and moves silver from 36.5% to a flattering 62.7%.

Report the mechanism grade; keep provenance as a column so a reader can see
which source supplied each ladder. Note that gold under the mechanism scheme
still rests on ChemAxon-derived ladders for most of its members — it means
well-founded protonation, not open provenance. Those are different claims.

**This is one axis.** A full evidence grade must also weigh the energy source
(measured anchor vs fitted prediction), the reported uncertainty, and whether
sources agree on direction. Presenting this as the whole grade would overstate
it.

## Licensing

`build_ladder_table.py` writes two tables. `literature_ladders.tsv` interleaves
Alberty with IUPAC values and is **gitignored**: the IUPAC Digitized pKa Dataset
is CC-BY-NC-4.0, which is incompatible with how this database is redistributed.
`literature_ladders_open.tsv` is Alberty only and is committed.

Citations are printable in both cases; it is the measured values that are
encumbered. The script expands each IUPAC row's reference code into its full
Serjeant/Perrin bibliographic entry, so those values are individually citable
even though they are not redistributed.

## Outputs

Written to `Biochemistry/Thermodynamics/ProtonationEvidence/`:

| file | contents |
|---|---|
| `reaction_evidence_grades.tsv` | one row per scored reaction: both grades, the limiting compound and why |
| `structural_match_classification.tsv` | per compound: safe match or refused, with the verdict |
| `ladder_requirements.tsv` | per head compound: the ladder shape the magnesium guard demands |
| `literature_ladders_open.tsv` | cited macroscopic ladders, Alberty only |

## The magnesium guard, stated usefully

The guard in `build_modelseed_cache.py` is written as a set containment. It
reduces exactly: with `n_h` the cache `atom_bag` proton count and `mg` the
proton counts the magnesium rows reference, a ladder is admissible iff it has at
least `max(0, n_h - min(mg))` values above pH 7 and at least
`max(0, max(mg) - n_h)` at or below it. Only the counts either side of pH 7
matter, not the values.

For several compounds — NAD, AMP, UDP, GTP — that demands more high-pKa steps
than the chemistry plausibly offers, and their magnesium rows skip proton counts
(UDP skips 10, GTP skips 11). A table indexed off one consistent ladder should
not have gaps, so those are likely defects in the magnesium data rather than
missing literature.
