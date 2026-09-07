#!/usr/bin/env python
"""The graphical abstract. NAR requires one; it is not optional.

SPEC, verbatim from latex/NAR_REQUIREMENTS.md section 4:
  5:2 aspect, 127x50mm / 5x2in MINIMUM, landscape, TIF/EPS/editable PDF,
  300-600 dpi, sans-serif 12-16pt, "be simple", "use text sparingly, mainly
  for labels", "read from top down or left to right", and it must "be original
  i.e. not an existing main or supplementary figure".

That last rule is why this is not a crop of figures 1-3. It is drawn from the
same measured numbers but says one thing they do not: the whole argument of the
paper in one left-to-right reading.

THE ARGUMENT: a curated network feeds predictors that DISAGREE -- on rxn00001
they disagree in sign -- so the release publishes every source with its own
uncertainty and grades the evidence, instead of promoting one number and hiding
the rest.

Numbers are IMPORTED from make_figures, not retyped. Transcribing them a second
time is how a graphical abstract ends up contradicting its own paper.

DELIBERATELY ABSENT: the gold-vs-bronze accuracy split. That claim is under a
circularity check (task #55) and is withheld from the abstract; it has no
business being the headline of the graphical abstract while that is unresolved.
The grade panel therefore shows how many reactions carry each grade, which is a
count and not in doubt.

Output: figures/graphical_abstract.pdf -- submitted as a SEPARATE file, never
\\includegraphics'd into main.tex.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
from pathlib import Path

import grace_style as grace
from make_figures import NUMBERS, BLUE, ORANGE, AQUA, VIOLET, NEUTRAL, k

OUT = Path(__file__).resolve().parent.parent / "figures" / "graphical_abstract.pdf"

# 7.5 x 3.0in = 190 x 76mm, comfortably over the 127x50mm floor, exactly 5:2.
W, H = 7.5, 3.0
GOLD_C, SILVER_C, BRONZE_C = "#b8860b", "#7d7d7d", "#a0522d"

# rxn00001, diphosphate phosphohydrolase -- the sign disagreement, in kcal/mol.
# Same values quoted in Scripts/Thermodynamics/reaction_direction_prompts.md.
RXN = [("Group contribution", 4.18, 2.24, AQUA),
       ("dGPredictor", -3.78, 1.48, VIOLET),
       ("eQuilibrator", -4.06, 0.05, BLUE)]


def main():
    plt.rcParams.update(grace.RC)
    fig = plt.figure(figsize=(W, H))
    fig.patch.set_facecolor("white")

    # LAYOUT. NAR mandates 12-16pt, which is very large relative to a 7.5in
    # canvas -- a first attempt with hand-placed labels overflowed the page and
    # overlapped all three panels. Two rules keep it honest:
    #   1. series names are AXIS TICK LABELS, so matplotlib places them and they
    #      cannot collide with a neighbouring panel;
    #   2. everything else sits in one of four disjoint horizontal bands.
    HEAD, SUB, AX_B, AX_H = 0.94, 0.82, 0.285, 0.375  # bands, figure fraction
    FOOT2 = 0.035

    cpd, rxn = NUMBERS["struct_total"][1], NUMBERS["energy_total"]
    have = NUMBERS["struct_total"][0]

    # ---------------- 1. the network ---------------------------------------
    fig.text(0.030, HEAD, "Curated network", fontsize=13, fontweight="bold",
             color="black", va="top")
    for i, (val, lab) in enumerate([(f"{cpd:,}", "compounds"),
                                    (f"{rxn:,}", "reactions"),
                                    (f"{have:,}", "structures")]):
        y = 0.635 - i * 0.135
        fig.text(0.030, y, val, fontsize=14, fontweight="bold", color=BLUE, va="center")
        fig.text(0.128, y, lab, fontsize=12, color="black", va="center")

    # ---------------- 2. the disagreement -----------------------------------
    ax = fig.add_axes([0.415, AX_B, 0.205, AX_H])
    for i, (name, v, e, col) in enumerate(RXN):
        ax.errorbar(v, len(RXN) - 1 - i, xerr=e, fmt="o", ms=7, color=col,
                    ecolor=col, elinewidth=2.0, capsize=3.5, capthick=2.0, zorder=4)
    ax.axvline(0, color="black", lw=1.1, zorder=2)
    ax.set_xlim(-7.5, 7.5); ax.set_ylim(-0.6, len(RXN) - 0.4)
    ax.set_yticks(range(len(RXN)))
    ax.set_yticklabels([n.replace(" ", "\n") for n, _, _, _ in RXN][::-1],
                       fontsize=12)
    for lbl, (name, _v, _e, col) in zip(ax.get_yticklabels()[::-1], RXN):
        lbl.set_color(col); lbl.set_fontweight("bold")
    ax.set_xticks([-5, 0, 5]); ax.tick_params(axis="x", labelsize=12)
    pass
    # Plain Unicode, not mathtext: NAR asks for a sans-serif figure, and
    # matplotlib renders $\\circ$ from Cmsy10, a SERIF Computer Modern symbol
    # font, which gets embedded alongside the sans text.
    ax.set_xlabel("\u0394G\u2032\u00b0  (kcal/mol)", fontsize=12, labelpad=0)
    grace.frame(ax)
    ax.tick_params(axis="y", length=0, pad=7)   # after frame(); see note above
    fig.text(0.290, HEAD, "Predictors disagree", fontsize=13, fontweight="bold",
             color="black", va="top")
    fig.text(0.290, SUB, "one reaction, three answers", fontsize=12,
             color="black", va="top")

    # ---------------- 3. graded, not promoted --------------------------------
    ax2 = fig.add_axes([0.845, AX_B, 0.128, AX_H])
    grades = [("GOLD", sum(NUMBERS["dir_grade"][0][1:]), GOLD_C),
              ("SILVER", sum(NUMBERS["dir_grade"][1][1:]), SILVER_C),
              ("BRONZE", sum(NUMBERS["dir_grade"][2][1:]), BRONZE_C)]
    for i, (lab, n, col) in enumerate(grades):
        ax2.barh(len(grades) - 1 - i, n, height=0.66, color=col,
                 edgecolor="black", lw=0.6, zorder=3)
        ax2.text(n + 700, len(grades) - 1 - i, k(n), va="center", fontsize=12,
                 color="black")
    ax2.set_xlim(0, 33000); ax2.set_ylim(-0.6, len(grades) - 0.4)
    ax2.set_yticks(range(len(grades)))
    ax2.set_yticklabels([g[0] for g in grades][::-1], fontsize=12,
                        fontweight="bold")
    for lbl, g in zip(ax2.get_yticklabels()[::-1], grades):
        lbl.set_color(g[2])
    ax2.set_xticks([])
    grace.frame(ax2, minor=None)
    ax2.tick_params(axis="y", length=0, pad=7)  # after frame(); see note above
    fig.text(0.690, HEAD, "Graded, not promoted", fontsize=13, fontweight="bold",
             color="black", va="top")
    fig.text(0.690, SUB, "every source published,", fontsize=12, color="black", va="top")
    fig.text(0.690, SUB - 0.085, "with its uncertainty", fontsize=12,
             color="black", va="top")

    # ---------------- flow --------------------------------------------------
    for x0 in (0.243, 0.638):
        fig.add_artist(FancyArrow(x0, AX_B + AX_H / 2, 0.028, 0, width=0.010,
                                  head_width=0.050, head_length=0.013,
                                  length_includes_head=True,
                                  color="#9a9a9a", transform=fig.transFigure))

    fig.text(0.030, FOOT2, "ModelSEED Biochemistry Database", fontsize=13,
             fontweight="bold", color="black", va="bottom")
    fig.text(0.472, FOOT2, "modelseed.org", fontsize=12, color=ORANGE, va="bottom")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, format="pdf")
    plt.close(fig)
    print(f"wrote {OUT}  ({W}x{H}in = {W*25.4:.0f}x{H*25.4:.0f}mm, {W/H:.2f}:1)")


if __name__ == "__main__":
    main()
