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
| `calibrate_sigma.py` | the empirical scale of each source's reported uncertainty |
| `figure_common.py` | shared palette, derived `NUMBERS` and panel helpers |
| `make_figure1_growth.py` | Figure 1, `../figures/figure1_growth.pdf` |
| `make_figureS1_pathway_classes.py` | Supplementary Figure S1, `../figures/figureS1_pathway_classes.pdf` — the pathway-class panel that was Figure 1D |
| `make_figure2_thermodynamics.py` | Figure 2, `../figures/figure2_thermodynamics.pdf` |
| `make_figure3_direction.py` | Figure 3, `../figures/figure3_direction.pdf` |
| `grace_style.py` | shared Grace/xmgrace visual style used by both |
| `pathway_distribution_of_growth.py` | where the 2020→2026 growth landed at MetaCyc class level — reviewer 1 comment 1 |

The graphical abstract is submitted as a **separate file** and is deliberately
never `\includegraphics`'d into `main.tex` — NAR requires it uploaded on its
own, and it must not duplicate a main figure. The generator script is not
tracked here; only the output (`../figures/graphical_abstract.{pdf,svg}`) is.

Findings written up from these live in `../data/`, dated -- untracked, local to the author's tree rather than in the repository.

## Running them

`calibrate_sigma.py` and `pmg_sensitivity.py` read caches and fitted parameters
from the eQuilibrator working tree, which is far too large to commit —
`data/` there is several gigabytes. That location is `$EQUILIBRATOR_DIR`,
defaulting to the analysis host's path. Paths into this repository derive from
the file location and need no configuration.

    EQUILIBRATOR_DIR=/path/to/eQuilibrator python3 pmg_sensitivity.py

They were moved here from that tree, where `ROOT` had been
`__file__/../..`. That is why `ROOT` is now named explicitly rather than
derived: the scripts sit in a different repository from the data they read, and
the relationship should be visible rather than implied by directory depth.

Every other script here reads only `Biochemistry/` and runs with no
`EQUILIBRATOR_DIR`.

## A caveat on `pmg_sensitivity`

This justifies a **live pipeline setting** — the release reports at pMg 3.0, and
that choice rests on this result. Filing the evidence under `Papers/` is right
by the rule above, but it puts the justification for a current setting in a
directory that plausibly gets archived after publication. If this condition
is ever revisited, start here.

## Removed

`crossval_pmg.py`, `evaluate_path_b.py`, `review_lost_reactions.py`,
`review_transport_and_llm.py` and `population_basis_table.py` were
investigation scripts whose output never made it into the manuscript (or, for
`review_transport_and_llm.py`, only half did — the LLM-ensemble numbers are in
S04, but the transport-cost breakdown was never cited). Removed rather than
kept as background, per the rule above: a script with no claim in the paper
isn't an analysis supporting the manuscript.
