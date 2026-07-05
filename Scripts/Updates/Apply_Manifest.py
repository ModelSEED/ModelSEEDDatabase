#!/usr/bin/env python
"""
Apply a single update manifest to the ModelSEEDDatabase.

Reads a YAML manifest from Biochemistry/Updates/, dispatches to the
matching handler, edits the appropriate files (or prints a dry-run
plan), and prints the cascade of pipeline stages that need to re-run.

Usage:
    python Apply_Manifest.py <manifest.yaml> [--dry-run] [--cascade]

Without --cascade, prints the cascade as a list of `Refresh_Pipeline.py`
arguments. With --cascade, runs Refresh_Pipeline.py automatically.

After successful (non-dry-run) application, the manifest is moved to
Biochemistry/Updates/applied/ with a YYYY-MM-DD prefix.

Exit codes:
  0   manifest applied (or dry-run plan printed)
  1   manifest invalid or apply failed
  2   environment error (PyYAML missing, file not found, etc.)
"""
import argparse
import datetime
import os
import shutil
import subprocess
import sys

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
BIOCHEM_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', 'Biochemistry'))
UPDATES_DIR  = os.path.join(BIOCHEM_ROOT, 'Updates')
APPLIED_DIR  = os.path.join(UPDATES_DIR, 'applied')

# Local imports
sys.path.insert(0, SCRIPT_DIR)
from cascade import stages_for, STAGES                 # noqa: E402
from handlers import HANDLERS                          # noqa: E402


REQUIRED_FIELDS = ('type', 'title', 'author', 'date', 'reason', 'target', 'change')


def die(msg, code=1):
    print(f'ERROR: {msg}', file=sys.stderr)
    sys.exit(code)


def load_manifest(path):
    try:
        import yaml
    except ImportError:
        die('PyYAML is required: pip install pyyaml', code=2)
    if not os.path.isfile(path):
        die(f'manifest not found: {path}', code=2)
    with open(path) as fh:
        m = yaml.safe_load(fh)
    if not isinstance(m, dict):
        die(f'manifest is not a YAML mapping: {path}')
    for f in REQUIRED_FIELDS:
        if f not in m:
            die(f'manifest missing required field: {f!r}')
    if m['type'] not in HANDLERS:
        die(f'unknown manifest type: {m["type"]!r} (handlers: {sorted(HANDLERS)})')
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('manifest', help='path to manifest YAML')
    ap.add_argument('--dry-run', action='store_true',
                    help='print what would change without writing files')
    ap.add_argument('--cascade', action='store_true',
                    help='also run Refresh_Pipeline.py for the resulting stages')
    args = ap.parse_args()

    m = load_manifest(args.manifest)
    print(f'manifest:    {args.manifest}')
    print(f'  type:      {m["type"]}')
    print(f'  title:     {m["title"]}')
    print(f'  author:    {m["author"]}  date: {m["date"]}')
    print()

    handler = HANDLERS[m['type']]
    ctx = {'biochem': BIOCHEM_ROOT, 'dry_run': args.dry_run}
    try:
        handler(m, ctx)
    except Exception as e:
        die(f'apply failed: {e}')

    cascade = stages_for(m['type'])
    print()
    print('Cascade (run in order):')
    for stage in cascade:
        label, _cmd = STAGES.get(stage, (f'<{stage}>', None))
        print(f'  -> {label}')

    if args.dry_run:
        print('\n(dry-run — no files were written, manifest not moved)')
        return

    if args.cascade:
        print()
        cmd = ['python3', os.path.join(SCRIPT_DIR, 'Refresh_Pipeline.py'),
               '--type', m['type']]
        print('Running:', ' '.join(cmd))
        rc = subprocess.call(cmd)
        if rc != 0:
            die(f'cascade failed (exit {rc})')

    # Move the manifest into applied/
    os.makedirs(APPLIED_DIR, exist_ok=True)
    today    = datetime.date.today().isoformat()
    base     = os.path.basename(args.manifest)
    if not base.startswith(today):
        base = f'{today}-{base}'
    target   = os.path.join(APPLIED_DIR, base)
    print(f'\narchive: {args.manifest} -> {target}')
    shutil.move(args.manifest, target)


if __name__ == '__main__':
    main()
