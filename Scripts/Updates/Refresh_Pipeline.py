#!/usr/bin/env python
"""
Refresh the pipeline stages for a given manifest type.

Apply_Manifest already knows which stages a given manifest type
triggers (see cascade.py). This script just runs them in order.

Usage:
    python Refresh_Pipeline.py --type <manifest_type>
    python Refresh_Pipeline.py --stages Print List UpdateStructures Reprint
    python Refresh_Pipeline.py --type structure_update --stop-after Reprint

Exit codes:
  0   all stages completed
  1   one or more stages failed (the script exits at the first failure)
"""
import argparse
import os
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from cascade import stages_for, STAGES                 # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument('--type',   help='manifest type; expand its cascade from cascade.py')
    grp.add_argument('--stages', nargs='+', help='explicit stage list (cascade keys)')
    ap.add_argument('--stop-after', help='stop after this stage (inclusive)')
    ap.add_argument('--dry-run', action='store_true',
                    help='print the planned stages and commands, run nothing')
    args = ap.parse_args()

    if args.type:
        try:
            stages = stages_for(args.type)
        except KeyError as e:
            print(f'ERROR: {e}', file=sys.stderr)
            sys.exit(1)
    else:
        stages = args.stages

    if args.stop_after:
        if args.stop_after not in stages:
            print(f'ERROR: --stop-after stage {args.stop_after!r} not in cascade', file=sys.stderr)
            sys.exit(1)
        stages = stages[:stages.index(args.stop_after) + 1]

    print('Plan:')
    for s in stages:
        label, cmd = STAGES[s]
        print(f'  {s}  ({label})  ->  {" ".join(cmd)}')

    if args.dry_run:
        return

    for s in stages:
        label, cmd = STAGES[s]
        print(f'\n--- {label} ---')
        rc = subprocess.call(cmd)
        if rc != 0:
            print(f'\nFAIL: stage {label!r} exited with {rc}', file=sys.stderr)
            sys.exit(1)
    print('\nAll stages completed.')


if __name__ == '__main__':
    main()
