#!/usr/bin/env python
"""Add measured pKa ladders to each compound's ``pkas`` dict as ``Literature``.

Added, never substituted. Marvin and MolGpKa entries stay exactly as they are,
so a user can see all three side by side and judge for themselves -- the same
policy the ``thermodynamics`` field follows. What the entries mean differs, and
that is the point of the ``kind`` flag:

* ``microscopic`` -- Marvin and MolGpKa: one value per ionizable site, computed
  on a single protonation state. Not a ladder. Phosphate reads 2.11 three times
  because its three oxygens are equivalent and a per-site predictor cannot see
  the charge left behind by the previous deprotonation.
* ``macroscopic`` -- Literature: successive dissociations of the whole
  molecule. Phosphate reads 7.22. This is what eQuilibrator's transform
  actually consumes.

Source is Alberty, *Thermodynamics of Biochemical Reactions* (Wiley 2003,
doi:10.1002/0471332607), via the ``BasicBiochemData3`` Mathematica package.
Values are derived from his species tables as
``pKa = [dfG(deprotonated) - dfG(protonated)] / (RT ln10)`` -- his own
``calcpK`` -- and validated against literature: ATP 7.60, Pi 7.22, PPi 9.46,
acetate 4.75, ammonia 9.25, citrate 6.39, succinate 5.64/4.21.

**The IUPAC dataset is deliberately not written here.** It reaches 1,453
compounds against Alberty's 31 and the pipeline does consume it, but it is
licensed CC-BY-NC-4.0 and this database is redistributed without that
restriction. It therefore stays in a gitignored flat file and never enters a
committed record. A consequence worth stating in any methods section: the
shipped energies were computed with IUPAC values in the resolution order, while
the ``Literature`` entries visible in these JSON files are Alberty only.

Two caveats on the values themselves. They are at **zero ionic strength**,
which is what eQuilibrator wants -- it applies Debye-Huckel itself from the
microspecies charge -- so no conversion is applied. And the ladders are
truncated to the physiologically interesting range, because Alberty tabulates
only species that matter between roughly pH 5 and 9: phosphate yields 7.22
alone, without its 2.15 and 12.35.
"""
import argparse
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# The eQuilibrator working tree, which supplies ``tools.modelseed_pkas``. It is
# not part of this repository and its location differs per host, so it is an
# environment variable rather than a literal. Previously hardcoded to one
# author's scratch directory, which made this script unrunnable for anyone else.
EQUILIBRATOR_DIR = os.environ.get("EQUILIBRATOR_DIR")
if EQUILIBRATOR_DIR:
    sys.path.insert(0, EQUILIBRATOR_DIR)

from pka_encoding import MACROSCOPIC, encode          # noqa: E402

BIOCHEM = HERE.parents[1] / "Biochemistry"
LABEL = "Literature"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--package", type=Path,
                    default=Path(os.environ.get(
                        "ALBERTY_PACKAGE", "BasicBiochemData3.m")),
                    help="Alberty's BasicBiochemData3.m. Not redistributable "
                         "with this repository; set ALBERTY_PACKAGE or pass "
                         "--package.")
    a = ap.parse_args()

    try:
        from tools import modelseed_pkas as msp
    except ImportError as exc:                       # fail with the reason, not a stack trace
        raise SystemExit(
            "cannot import tools.modelseed_pkas -- set EQUILIBRATOR_DIR to the "
            f"eQuilibrator working tree (currently {EQUILIBRATOR_DIR or 'unset'}). "
            f"Original error: {exc}")

    ladders = msp.load_alberty_pkas(a.package)
    if not ladders:
        raise SystemExit(f"no ladders derived from {a.package}")
    print(f"Alberty ladders: {len(ladders)}")

    seen = {"written": 0, "unchanged": 0, "removed": 0, "no such compound": 0}
    found = set()
    for path in sorted(BIOCHEM.glob("compound_*.json")):
        entries = json.load(path.open())
        touched = False
        for e in entries:
            pkas = e.get("pkas")
            if e["id"] in ladders:
                found.add(e["id"])
                rec = {"kind": MACROSCOPIC,
                       "pKa": encode(ladders[e["id"]]),
                       "pKb": ""}
                pkas = e.setdefault("pkas", {})
                if pkas.get(LABEL) == rec:
                    seen["unchanged"] += 1
                else:
                    pkas[LABEL] = rec
                    seen["written"] += 1
                    touched = True
            elif isinstance(pkas, dict) and LABEL in pkas:
                # a compound dropped from the source must lose the entry too,
                # or the record keeps asserting a value nothing stands behind
                del pkas[LABEL]
                seen["removed"] += 1
                touched = True
        if touched and not a.dry_run:
            with path.open("w") as fh:
                fh.write(json.dumps(entries, indent=4, sort_keys=True))

    seen["no such compound"] = len(set(ladders) - found)
    for k, v in seen.items():
        print(f"  {k:<20}{v:>7,}")
    if seen["no such compound"]:
        print(f"    unmapped: {sorted(set(ladders) - found)}")
    if a.dry_run:
        print("\n(dry run - nothing written)")


if __name__ == "__main__":
    main()
