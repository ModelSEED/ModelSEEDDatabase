#!/usr/bin/env python
"""Grace/xmgrace visual style for matplotlib.

The 2020 ModelSEED paper's figures were drawn with PyGrace, and the look Sam
wants back is Grace's, not PyGrace's API. That look is structural, not chromatic,
so it reproduces exactly in matplotlib:

  * a closed black frame on all four sides, not two open spines
  * major AND minor ticks pointing INWARD, on all four sides
  * no grid -- Grace ships with the grid off, and the frame plus inward ticks
    carry the reading
  * Helvetica throughout
  * pure white canvas, black rules, labels outside the frame

WHAT IS DELIBERATELY NOT COPIED: Grace's default colour cycle
(black/red/green/blue/yellow...). Its red and green are indistinguishable to a
deuteranopic reader, and the palette these figures already use was validated for
CVD separation. Geometry and typography are Grace's; the hues stay accessible.

Helvetica comes from TeX Gyre Heros, the URW Nimbus Sans clone shipped with TeX
Live -- the same metrics Grace used. If TeX Live moves, FONT falls back to
DejaVu Sans and everything still renders, just not in Helvetica.
"""
from pathlib import Path

import matplotlib.font_manager as fm

TEXLIVE_HEROS = Path("/scratch/seaver/texlive/2026/texmf-dist/fonts/opentype"
                     "/public/tex-gyre")

FRAME = "#000000"   # Grace draws the frame and ticks in pure black
INK = "#000000"


def _register_helvetica():
    """Return the family name to use, registering TeX Gyre Heros if present."""
    if not TEXLIVE_HEROS.is_dir():
        return "DejaVu Sans"
    added = 0
    for otf in sorted(TEXLIVE_HEROS.glob("texgyreheros-*.otf")):
        try:
            fm.fontManager.addfont(str(otf))
            added += 1
        except Exception:
            pass
    if not added:
        return "DejaVu Sans"
    names = {f.name for f in fm.fontManager.ttflist}
    return "TeX Gyre Heros" if "TeX Gyre Heros" in names else "DejaVu Sans"


FONT = _register_helvetica()

#: rcParams implementing the Grace look. Apply with plt.rcParams.update(RC).
RC = {
    # A list, not a string: matplotlib falls back PER GLYPH, and TeX Gyre Heros
    # has no arrows (U+2190/2192/2194), which figure 3 uses as direction marks.
    # Helvetica carries the text; DejaVu supplies only the glyphs Heros lacks.
    "font.family": [FONT, "DejaVu Sans"],
    "font.size": 7.2,
    "text.color": INK,

    # Closed black frame. Grace's default frame is a box, and its weight is
    # visually heavier than matplotlib's default hairline.
    "axes.edgecolor": FRAME,
    "axes.linewidth": 1.0,
    "axes.labelcolor": INK,
    "axes.grid": False,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "axes.spines.left": True,
    "axes.spines.bottom": True,

    # Inward ticks on all four sides -- the single most recognisable Grace cue.
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 5.0, "ytick.major.size": 5.0,
    "xtick.minor.size": 2.6, "ytick.minor.size": 2.6,
    "xtick.major.width": 0.9, "ytick.major.width": 0.9,
    "xtick.minor.width": 0.7, "ytick.minor.width": 0.7,
    "xtick.color": FRAME, "ytick.color": FRAME,
    "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "xtick.labelsize": 6.8, "ytick.labelsize": 6.8,
    "xtick.major.pad": 2.5, "ytick.major.pad": 2.5,

    # Grace legends sit in a square white box with a thin black border.
    "legend.frameon": True,
    "legend.edgecolor": FRAME,
    "legend.facecolor": "white",
    "legend.fancybox": False,
    "legend.framealpha": 1.0,
    "legend.borderpad": 0.4,
    "legend.fontsize": 6.4,

    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,   # embed as TrueType so the PDF is editable downstream
}


def frame(ax, minor="value", value_axis="x"):
    """Apply the Grace frame to one axes.

    minor: "value" adds minor ticks to the measured axis only, "both", or None.
           Category axes get none -- a minor tick between two bar rows means
           nothing.
    """
    from matplotlib.ticker import AutoMinorLocator

    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(True)
        ax.spines[s].set_color(FRAME)
        ax.spines[s].set_linewidth(1.0)

    ax.tick_params(which="major", direction="in", length=5.0, width=0.9,
                   color=FRAME, top=True, right=True, labelcolor=INK)
    ax.tick_params(which="minor", direction="in", length=2.6, width=0.7,
                   color=FRAME, top=True, right=True)

    if minor in ("value", "both"):
        if value_axis == "x" or minor == "both":
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        if value_axis == "y" or minor == "both":
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.grid(False)
    ax.set_axisbelow(True)
    return ax
