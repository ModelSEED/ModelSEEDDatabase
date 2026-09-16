#!/usr/bin/env python
"""Convert the hand-authored graphical abstract from SVG to submission PDF.

NAR requires the graphical abstract as TIF, EPS or *editable* PDF -- SVG is not
accepted -- so the shipped artefact is a conversion, not a drawing. The artwork
itself is hand-authored in ../figures/graphical_abstract.svg; edit that, then
re-run this. Nothing here draws anything.

The spec, from the Database Issue guidelines (see latex/NAR_REQUIREMENTS.md §4):
5:2 aspect ratio, at least 127x50 mm, landscape, 300-600 dpi, sans-serif font.
The SVG is authored at 720x288 pt (10x4 in, 254x101.6 mm), which is exactly 5:2
and twice the minimum, so the conversion preserves compliance rather than
establishing it. This script asserts the geometry so a later edit to the SVG
cannot silently break it.

Requires cairosvg (pip install cairosvg; needs libcairo present).
"""
import sys
from pathlib import Path

FIGDIR = Path(__file__).resolve().parent.parent / "figures"
SRC = FIGDIR / "graphical_abstract.svg"
OUT = FIGDIR / "graphical_abstract.pdf"


def main():
    try:
        import cairosvg
    except ImportError:
        sys.exit("cairosvg is not installed: pip install cairosvg")
    if not SRC.exists():
        sys.exit(f"missing source artwork: {SRC}")
    cairosvg.svg2pdf(url=str(SRC), write_to=str(OUT))

    # Geometry is a submission requirement, so verify rather than assume.
    try:
        from pypdf import PdfReader
    except ImportError:
        print(f"wrote {OUT} (install pypdf to verify geometry)")
        return
    box = PdfReader(str(OUT)).pages[0].mediabox
    w, h = float(box.width), float(box.height)
    mm_w, mm_h = w / 72 * 25.4, h / 72 * 25.4
    assert abs(w / h - 2.5) < 0.01, f"aspect {w/h:.3f}, NAR requires 5:2"
    assert mm_w >= 127 and mm_h >= 50, f"{mm_w:.0f}x{mm_h:.0f} mm below the 127x50 minimum"
    assert w > h, "must be landscape"
    print(f"wrote {OUT}  ({mm_w:.0f}x{mm_h:.0f} mm, aspect {w/h:.3f})")


if __name__ == "__main__":
    main()
