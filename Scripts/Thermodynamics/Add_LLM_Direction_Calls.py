#!/usr/bin/env python

if __name__ == "__main__":
    # Validate arguments BEFORE importing anything or touching the database.
    import argparse as _argparse
    _ap = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter)
    _ap.add_argument("--calls", default="/scratch/jplfaria/share/"
                     "reaction_direction_council/variant_A_directions.tsv",
                     help="council output TSV: reaction_id, council_direction")
    _ap.add_argument("--skip-abstentions", action="store_true",
                     help="omit reactions the council declined ('?') instead of "
                          "recording the declination")
    _ap.add_argument("--remove", action="store_true",
                     help="delete every LLMs entry and write nothing")
    _ap.add_argument("--dry-run", action="store_true")
    _ARGS = _ap.parse_args()

#
# Records the LLM ensemble's direction call as a fourth thermodynamics source,
# under the label "LLMs".
#
# The entry deliberately carries NO energy:
#
#     "LLMs": ["", "", ">"]
#
# The ensemble reasons from the reaction rather than from a number, so there is
# no dG and no uncertainty to report -- only the direction. Blank is the honest
# encoding, and it is also what keeps the call out of everything that consumes
# energies:
#
#   * SourceGrading/optimize_thermo_source_assignment.py admits a source only
#     when float(val[0]) parses, so "" excludes LLMs from grading. This is the
#     stated policy: the calls ship in the database but do not enter the
#     evidence grades.
#   * Estimate_Reaction_Reversibility.py iterates DB_LEVELS = (EQ, GC, DGP), so
#     the cascade never sees LLMs.
#   * Apply_Evidence_Grades_And_Recommendation.py has a fixed PRECEDENCE of
#     eQuilibrator > dGPredictor > Group contribution, so LLMs cannot become
#     the recommended direction in the `reversibility` field.
#   * Compile_Biochemistry_for_SOLR.py wraps float() in try/except, so the
#     child doc indexes with `operator` and no `energy`/`error`.
#
# Two SOLR-side effects DO follow, because both are computed over the whole
# thermodynamics dict: n_sources_thermodynamics gains one for every reaction
# with a call, and sources_agree_direction starts taking the LLM call into
# account. Neither is wrong, but both change meaning, so they are reported.
#
# Idempotent: re-running overwrites the LLMs entry and touches nothing else.
#
import csv
import os
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '..', '..', 'Libs', 'Python'))
from BiochemPy import Reactions

LABEL = "LLMs"
VALID = ('>', '<', '=', '?')


def main(args):
    reactions_helper = Reactions()
    reactions_dict = reactions_helper.loadReactions()

    if args.remove:
        n = 0
        for rxn in reactions_dict.values():
            t = rxn.get('thermodynamics')
            if isinstance(t, dict) and LABEL in t:
                del t[LABEL]
                n += 1
        print(f"removed {LABEL} from {n:,} reactions")
        if not args.dry_run:
            reactions_helper.saveReactions(reactions_dict)
        return 0

    if not os.path.exists(args.calls):
        print(f"ERROR: {args.calls} not found", file=sys.stderr)
        return 1

    calls, malformed = {}, 0
    with open(args.calls) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            rid = (row.get('reaction_id') or '').strip()
            d = (row.get('council_direction') or '').strip()
            if not rid or d not in VALID:
                malformed += 1
                continue
            calls[rid] = d
    print(f"read {len(calls):,} calls from {os.path.basename(args.calls)}"
          f"  ({malformed} malformed)")

    unknown = sorted(set(calls) - set(reactions_dict))
    if unknown:
        print(f"WARNING: {len(unknown):,} reaction ids are not in the database, "
              f"skipped (e.g. {unknown[:3]})")

    written = skipped = replaced = 0
    by_dir = {d: 0 for d in VALID}
    for rid, d in calls.items():
        rxn = reactions_dict.get(rid)
        if rxn is None:
            continue
        if d == '?' and args.skip_abstentions:
            skipped += 1
            continue
        t = rxn.get('thermodynamics')
        if not isinstance(t, dict):
            t = rxn['thermodynamics'] = {}
        if LABEL in t:
            replaced += 1
        t[LABEL] = ["", "", d]
        written += 1
        by_dir[d] += 1

    print(f"wrote {LABEL} on {written:,} reactions "
          f"({replaced:,} replaced an existing entry, {skipped:,} abstentions skipped)")
    print("  by direction: " + ", ".join(f"{d} {by_dir[d]:,}" for d in VALID))
    covered = sum(1 for r in reactions_dict.values()
                  if isinstance(r.get('thermodynamics'), dict)
                  and LABEL in r['thermodynamics'])
    print(f"  reactions now carrying an {LABEL} entry: {covered:,} "
          f"of {len(reactions_dict):,}")

    if args.dry_run:
        print("dry run — nothing written")
        return 0
    print("Saving reactions")
    reactions_helper.saveReactions(reactions_dict)
    return 0


if __name__ == "__main__":
    sys.exit(main(_ARGS))
