#!/usr/bin/env python
"""Apply the 2020 reversibility policy: eQuilibrator primary, GC fallback.

WHY THIS EXISTS. Estimate_Reaction_Reversibility.py writes the TOP-LEVEL
`reversibility` field from whichever source you name, and it rewrites every
reaction -- there is no scoping option. So the sources are NOT additive: run two
in sequence and the second silently replaces the first everywhere. Running all
of them in order ends with DGPM, for which no data exists, and sets all 56,012
reactions to '?'. That is a whole-database wipe from a plausible-looking loop,
and it is why this policy belongs in a script rather than in someone's shell
history.

THE POLICY, reproduced from commits 636e16fa and bc944f90:

  1. eQuilibrator energies under the HISTORICAL GC RULE SET
     (EQ --heuristics GC, what the estimator calls "old behaviour").
     EQ_HEURISTICS is Beber 2022 over the Noor 2012 index and postdates the
     2020 paper, so it is deliberately not used.
  2. Group contribution as fallback for everything eQuilibrator cannot reach,
     applied by snapshotting the eQ calls, running GC wholesale, and restoring
     the eQ-covered reactions.

Only the inputs change between releases; the rules are fixed.
"""
import argparse
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
BIOCHEM = HERE.parents[1] / "Biochemistry"
EST = HERE / "Estimate_Reaction_Reversibility.py"


def read_state():
    """{rxn: (reversibility, has_eq_energy)} across every reaction shard."""
    out = {}
    for path in sorted(BIOCHEM.glob("reaction_*.json")):
        for r in json.load(path.open()):
            eq = (r.get("thermodynamics") or {}).get("eQuilibrator")
            out[r["id"]] = (r.get("reversibility"),
                            isinstance(eq, list) and eq and eq[0] is not None)
    return out


def write_reversibility(calls):
    """Restore `reversibility` for the given {rxn: value}, leaving others alone."""
    for path in sorted(BIOCHEM.glob("reaction_*.json")):
        entries = json.load(path.open())
        touched = False
        for r in entries:
            if r["id"] in calls and r.get("reversibility") != calls[r["id"]]:
                r["reversibility"] = calls[r["id"]]
                touched = True
        if touched:
            path.open("w").write(json.dumps(entries, indent=4, sort_keys=True))


def run(*args):
    p = subprocess.run([sys.executable, str(EST), *args], cwd=HERE,
                       capture_output=True, text=True)
    if p.returncode:
        sys.exit(f"{' '.join(args) or 'top'} failed:\n{p.stdout[-2000:]}{p.stderr[-2000:]}")
    return p.stdout


def census(state):
    return dict(Counter(v for v, _ in state.values()).most_common())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="report what the two stages would do; restores the "
                         "starting state afterwards")
    a = ap.parse_args()

    before = read_state()
    print(f"before         {census(before)}")

    run("EQ", "--heuristics", "GC")
    after_eq = read_state()
    eq_calls = {r: v for r, (v, has_eq) in after_eq.items() if has_eq}
    print(f"after eQ       {census(after_eq)}   ({len(eq_calls)} eQ-covered)")

    run("GC")
    print(f"after GC       {census(read_state())}")

    write_reversibility(eq_calls)
    final = read_state()
    print(f"after restore  {census(final)}")

    changed = [r for r in before if before[r][0] != final[r][0]]
    trans = Counter((before[r][0], final[r][0]) for r in changed)
    print(f"\nchanged: {len(changed)}")
    for (x, y), n in trans.most_common(12):
        print(f"  {x!s:>4} -> {y!s:<4} {n:6d}")

    # The restore must not have let GC overwrite an eQ-covered reaction.
    bad = [r for r, v in eq_calls.items() if final[r][0] != v]
    print(f"\neQ-covered reactions overwritten by the GC pass: {len(bad)}  "
          f"{'OK' if not bad else 'BUG'}")

    if a.dry_run:
        write_reversibility({r: v for r, (v, _) in before.items()})
        print("\n(dry run - starting state restored)")


if __name__ == "__main__":
    main()
