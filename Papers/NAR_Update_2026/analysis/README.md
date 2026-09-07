# Analyses supporting the manuscript

Scripts here exist to produce a number, table or argument in the paper. They are
**not** part of the pipeline that builds the database, and nothing under
`Biochemistry/` depends on them.

## The rule

> Does the script's output ship in `Biochemistry/`, or does it only ever appear
> in the manuscript?

Output destination, not subject matter. Most of the pipeline is "analysis" in
the ordinary sense, so that word does not discriminate; where the artefact lands
does, and it stays checkable as things move.

* ships in `Biochemistry/` → `Scripts/`, it is pipeline
* appears only in the paper → here

Applying it moved two things back that a looser reading would have taken:
`check_mobile_h.py` looks like a one-off investigation, but
`grade_protonation_evidence.py` consumes its
`structural_match_classification.tsv`, and `seed_mapping.tsv` cannot substitute
because it does not distinguish "no candidate found" from "candidate refused on
tautomer grounds". `ladder_requirements.tsv` is purely diagnostic, but splitting
one script's two outputs across two trees costs more than the tidiness is worth.
Both stayed in `Scripts/Thermodynamics/ProtonationEvidence/`.

## Contents

| script | supports |
|---|---|
| `pmg_sensitivity.py` | how far dG'° moves with the magnesium condition — the "Reported conditions" argument in the thermodynamics methods |
| `crossval_pmg.py` | whether pMg 14 is the right default for the 66% of TECRDB rows that never recorded magnesium; the training/prediction inconsistency table |
| `calibrate_sigma.py` | the empirical scale of each source's reported uncertainty |
| `evaluate_path_b.py` | cache variants compared over real ModelSEED reactions |
| `review_lost_reactions.py` | reactions upstream eQuilibrator can compute and our build cannot |
| `make_figures.py` | the three main-text figures, `../figures/main_figures_draft.pdf` |
| `make_graphical_abstract.py` | the mandatory graphical abstract, `../figures/graphical_abstract.pdf` |
| `grace_style.py` | shared Grace/xmgrace visual style used by both |

The graphical abstract is submitted as a **separate file** and is deliberately
never `\includegraphics`'d into `main.tex` — NAR requires it uploaded on its
own, and it must not duplicate a main figure.

Findings written up from these live in `../data/`, dated -- untracked, local to the author's tree rather than in the repository.

## Running them

They read caches and fitted parameters from the eQuilibrator working tree, which
is far too large to commit — `data/` there is several gigabytes. That location is
`$EQUILIBRATOR_DIR`, defaulting to the analysis host's path. Paths into this
repository derive from the file location and need no configuration.

    EQUILIBRATOR_DIR=/path/to/eQuilibrator python3 pmg_sensitivity.py

They were moved here from that tree, where `ROOT` had been
`__file__/../..`. That is why `ROOT` is now named explicitly rather than
derived: the scripts sit in a different repository from the data they read, and
the relationship should be visible rather than implied by directory depth.

## A caveat on `pmg_sensitivity` and `crossval_pmg`

These justify **live pipeline settings** — the release reports at pMg 3.0, and
that choice rests on these results. Filing the evidence under `Papers/` is right
by the rule above, but it puts the justification for a current setting in a
directory that plausibly gets archived after publication. If these conditions
are ever revisited, start here.
