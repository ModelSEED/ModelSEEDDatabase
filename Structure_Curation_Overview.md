# Resolving a structure conflict

Most compounds in this database have one candidate structure and need no
judgement: 35,566 of 37,013 picks are `single_structure`. This document is
about the remainder — compounds where the source databases disagree and
somebody has to choose.

You resolve a conflict by **adding one line to a file and opening a pull
request**. You do not run the pipeline, and you do not edit the compound
records. Both of those are our job during review.

---

## 1. Find a conflict

Two reports, both regenerated whenever structures change, both tab-separated
and **without a header row**:

| file | rows | compounds | columns |
|---|---:|---:|---|
| `Biochemistry/Structures/_reports/Structure_Conflicts.txt` | 2,999 | 1,028 | `cpd_id · type · stage · structure · source_id · source_db` |
| `Biochemistry/Structures/_reports/Formula_Conflicts.txt` | 143 | 56 | `cpd_id · type · stage · formula · charge · source_id · source_db` |

One line per candidate, so a compound with three disagreeing sources appears
three times. Structure conflicts are disagreements about the molecule;
formula conflicts are disagreements about its composition or charge, which
are usually the more consequential because they break mass balance.

`Biochemistry/Structures/Pick_Reasons.txt` says what the pipeline currently
does with each compound. The 70 lines reading `formula_conflict_no_pick` are
the cases it refuses outright — those are the most valuable to resolve,
because right now the compound has no structure at all.

Being listed as a conflict does not mean the pipeline made a bad choice. It
resolves most of them automatically, by source agreement, by preferring the
stereochemically richer structure, or by source priority
(MetaCyc → KEGG → ChEBI → Rhea). A manual pick is worth making when you know
something the ordering cannot: which tautomer is physiological, which charge
state belongs at pH 7, which source simply has it wrong.

## 2. Record your pick

Add a row to **your own file** at

```
Biochemistry/Curation/overrides/structure_picks/<yourname>.tsv
```

Create it if it does not exist; the filename is your curator identity and it
is published with the pick. Seven tab-separated columns, header included:

```
cpd_id	format	structure	source_db	source_id	date	rationale
```

```
cpd00074	InChI	InChI=1S/S	KEGG	Original	2026-07-03	Elemental sulfur; KEGG's neutral form matches the formula S/0 carried by the compound record
```

- **`format`** — `InChI` or `SMILE`
- **`structure`** — the literal string you are choosing, usually copied from
  one of the conflict rows
- **`source_db` / `source_id`** — where it came from, e.g. `KEGG` /
  `Original`. Use `curator` / `manual` if you are supplying a structure that
  is in none of the sources
- **`rationale`** — free text, and the part that matters most. Say *why*,
  not *what*: "beta anomer is the physiological form" is useful, "picked
  MetaCyc" is not. It ships in the released data

One row per compound. Do not edit another curator's file — if two files claim
the same compound, the pipeline warns and takes one arbitrarily.

## 3. Open a pull request

Fork the repository, commit your file on a branch, and open a PR against
`dev`. This assumes you are comfortable with git and GitHub; we have not
written that part down.

**Please do not run the update pipeline and do not commit regenerated
outputs.** A single pick can change thousands of derived lines across the
structure reports, the compound records and both energy layers, which makes
the PR unreviewable and buries the one line that actually matters. A PR that
adds `N` lines to one TSV is one we can read in a minute.

## 4. What happens next

You do not run any of this yourself, which is why no script is named anywhere in
this guide. Validation is a maintainer step: it happens on our side, against the
whole database, once you open the pull request. Running the cascade against a
partial checkout would not see every reaction the compound appears in, so its
mass-balance verdict would be misleading, and that is exactly the check a pick
most needs.

We fetch the branch, apply the pick, and run the cascade
(`UpdateStructures → Reprint → BuildProvenance`) on our side. Then we tell you
what it did. Typically:

- the compound's `formula`, `charge`, `inchikey` and `smiles`, all re-derived
  through RDKit from the structure you chose, so they stay internally
  consistent rather than being copied from the source;
- whether the change fixes or breaks mass balance in any reaction the
  compound appears in — the usual reason a pick gets sent back;
- the line that will appear in `Pick_Reasons.txt`, reading
  `manual_curation:<yourname>`.

Two things worth knowing about how the pipeline treats a pick. It
**short-circuits the entire automatic cascade** — no agreement test, no
stereochemistry comparison, no source priority. And the underlying
disagreement is **still written to the conflict report**, so an override
records a decision without concealing the evidence for it.

If RDKit cannot parse your structure, the pick is skipped with a warning and
the automatic choice stands. That failure is silent in the data, which is
another reason we test before merging rather than after.
