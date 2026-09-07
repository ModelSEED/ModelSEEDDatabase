"""
Per-type handlers for Apply_Manifest. Each handler is a function that
takes the parsed manifest dict and a context (BIOCHEM_ROOT path,
dry_run flag) and either edits the right files or prints what it
would change.

Handlers are intentionally small — they do the file edit, nothing
else. The cascade (which scripts to re-run after) lives in cascade.py
so a handler stays focused on "what is the minimal data change for
this manifest type?".

Adding a new handler: define `apply_<type>(manifest, ctx)`, register
it in HANDLERS at the bottom.
"""

if __name__ == "__main__":
    # Validate arguments BEFORE importing anything or touching the database.
    # These scripts mutate the database, and without this an unknown flag or a
    # mistyped mode was silently ignored and the script ran with its defaults:
    # asking Estimate_Reaction_Reversibility.py for --help rewrote 122 files.
    # Placed above the imports so --help works even where a dependency is
    # missing from the path.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import csv
import os
import shutil


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _read_tsv(path):
    if not os.path.isfile(path):
        return [], None
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith('#')]
    reader = csv.DictReader(lines, dialect='excel-tab')
    return list(reader), reader.fieldnames


def _write_tsv(path, fieldnames, rows):
    with open(path, 'w') as fh:
        fh.write('\t'.join(fieldnames) + '\n')
        for row in rows:
            fh.write('\t'.join(str(row.get(f, '') or '') for f in fieldnames) + '\n')


def _say(ctx, msg):
    prefix = '[DRY-RUN] ' if ctx['dry_run'] else '          '
    print(prefix + msg)


# ---------------------------------------------------------------------------
# structure_update: edit one row in inchi.tsv / smiles.tsv / inchikey.tsv
# ---------------------------------------------------------------------------

def apply_structure_update(manifest, ctx):
    src    = manifest['target']['source']
    ext_id = manifest['target']['external_id']
    field  = manifest['target']['field']
    new    = manifest['change']['new_value']

    if field not in ('inchi', 'smiles', 'inchikey'):
        raise ValueError(f"unsupported field {field!r} (use inchi/smiles/inchikey)")
    fname = f'{field}.tsv'
    path  = os.path.join(ctx['biochem'], 'Structures', src, fname)
    rows, headers = _read_tsv(path)
    if not headers:
        raise FileNotFoundError(f"missing {path}")

    matched = 0
    for row in rows:
        if row['external_id'] == ext_id:
            old = row.get(field, '')
            if old == new:
                _say(ctx, f"{src}/{fname}: {ext_id}.{field} already equals new value, no change")
            else:
                _say(ctx, f"{src}/{fname}: {ext_id}.{field}\n             old: {old}\n             new: {new}")
                row[field] = new
                # Clear formula/charge so Print refreshes them; if we
                # leave them stale, downstream consumers see a structure
                # that doesn't match its formula until Print runs.
                if field in ('inchi', 'smiles'):
                    row['formula'] = ''
                    row['charge']  = ''
            matched += 1

    if matched == 0:
        raise KeyError(f"{ext_id} not found in {path}; add a row first or check the source")

    if not ctx['dry_run']:
        _write_tsv(path, headers, rows)


# ---------------------------------------------------------------------------
# protonation_replace: drop a new <tool>_<ver>_ph<n>.tsv file in place
# ---------------------------------------------------------------------------

def apply_protonation_replace(manifest, ctx):
    src     = manifest['target']['source']
    src_path  = manifest['change']['source_file']
    dest_rel  = manifest['change']['destination']
    dest_path = os.path.join(ctx['biochem'], 'Structures', src, dest_rel)

    if not os.path.isfile(src_path):
        raise FileNotFoundError(f"source file {src_path} not found")

    _say(ctx, f"copy {src_path} -> {dest_path}")
    if not ctx['dry_run']:
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        shutil.copy(src_path, dest_path)

    for retire_rel in manifest['change'].get('retire_existing', []) or []:
        retire_path = os.path.join(ctx['biochem'], 'Structures', src, retire_rel)
        if os.path.isfile(retire_path):
            _say(ctx, f"retire (delete) {retire_path}")
            if not ctx['dry_run']:
                os.remove(retire_path)


# ---------------------------------------------------------------------------
# override_add: append a row to Curation/overrides/<file>
# ---------------------------------------------------------------------------

def apply_override_add(manifest, ctx):
    rel  = manifest['target']['file']
    path = os.path.join(ctx['biochem'], rel)
    rows, headers = _read_tsv(path)
    if headers is None:
        raise FileNotFoundError(f"override file {path} not found")

    new_row = manifest['change']['add_row']
    if any(r.get('ID') == new_row.get('ID') for r in rows):
        _say(ctx, f"{rel}: id {new_row.get('ID')} already present, no change")
        return

    _say(ctx, f"{rel}: append row {new_row}")
    if not ctx['dry_run']:
        rows.append({k: str(v) for k, v in new_row.items()})
        _write_tsv(path, headers, rows)


# ---------------------------------------------------------------------------
# ignore_add: append a row to Curation/ignores/<file>
# ---------------------------------------------------------------------------

def apply_ignore_add(manifest, ctx):
    rel    = manifest['target']['file']
    path   = os.path.join(ctx['biochem'], rel)
    create = manifest['target'].get('create_if_missing', False)
    new    = manifest['change']['add_row']

    if not os.path.isfile(path):
        if not create:
            raise FileNotFoundError(f"{path} not found and create_if_missing is false")
        _say(ctx, f"create new ignore file {path}")
        if not ctx['dry_run']:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, 'w') as fh:
                fh.write('ID\tAccepted\tNotes\n')

    _say(ctx, f"{rel}: append {new['external_id']}\t{new.get('accepted','None')}\t{new.get('notes','')}")
    if not ctx['dry_run']:
        with open(path, 'a') as fh:
            fh.write(f"{new['external_id']}\t{new.get('accepted','None')}\t{new.get('notes','')}\n")


# ---------------------------------------------------------------------------
# alias_add / alias_remove: edit Aliases/Unique_ModelSEED_Compound_Aliases.txt
# ---------------------------------------------------------------------------

def _alias_path(ctx):
    return os.path.join(ctx['biochem'], 'Aliases', 'Unique_ModelSEED_Compound_Aliases.txt')


def apply_alias_add(manifest, ctx):
    cpd  = manifest['target']['modelseed_id']
    path = _alias_path(ctx)
    rows, headers = _read_tsv(path)
    if headers is None:
        raise FileNotFoundError(path)

    # The alias file uses 'ModelSEED ID', 'External ID', 'Source' columns
    existing = {(r.get('ModelSEED ID'), r.get('External ID'), r.get('Source')) for r in rows}
    added = 0
    for entry in manifest['change']['add']:
        ext = str(entry['external_id'])
        src = entry['source']
        key = (cpd, ext, src)
        if key in existing:
            _say(ctx, f"alias already present: {cpd} -> {src}:{ext}")
            continue
        _say(ctx, f"alias add: {cpd} -> {src}:{ext}")
        rows.append({'ModelSEED ID': cpd, 'External ID': ext, 'Source': src})
        existing.add(key)
        added += 1
    if added and not ctx['dry_run']:
        _write_tsv(path, headers, rows)


def apply_alias_remove(manifest, ctx):
    cpd  = manifest['target']['modelseed_id']
    path = _alias_path(ctx)
    rows, headers = _read_tsv(path)
    if headers is None:
        raise FileNotFoundError(path)

    drop = {(str(e['external_id']), e['source']) for e in manifest['change']['remove']}
    keep = []
    removed = 0
    for r in rows:
        if r.get('ModelSEED ID') == cpd and (r.get('External ID'), r.get('Source')) in drop:
            _say(ctx, f"alias remove: {cpd} -> {r.get('Source')}:{r.get('External ID')}")
            removed += 1
            continue
        keep.append(r)
    if removed and not ctx['dry_run']:
        _write_tsv(path, headers, keep)


# ---------------------------------------------------------------------------
# pka_replace: replace one <source>/pkas/<file>.tsv
# ---------------------------------------------------------------------------

def apply_pka_replace(manifest, ctx):
    src       = manifest['target']['source']
    src_path  = manifest['change']['source_file']
    dest_rel  = manifest['change']['destination']
    dest_path = os.path.join(ctx['biochem'], 'Structures', src, dest_rel)
    if not os.path.isfile(src_path):
        raise FileNotFoundError(src_path)
    _say(ctx, f"copy {src_path} -> {dest_path}")
    if not ctx['dry_run']:
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        shutil.copy(src_path, dest_path)


# ---------------------------------------------------------------------------
# dispatch
# ---------------------------------------------------------------------------

HANDLERS = {
    'structure_update':    apply_structure_update,
    'protonation_replace': apply_protonation_replace,
    'override_add':        apply_override_add,
    'ignore_add':          apply_ignore_add,
    'alias_add':           apply_alias_add,
    'alias_remove':        apply_alias_remove,
    'pka_replace':         apply_pka_replace,
}
