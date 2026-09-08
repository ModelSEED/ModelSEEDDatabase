#!/usr/bin/env python
"""Regenerate the per-source Marvin pKa bundles with ChemAxon cxcalc.

Replaces the pKaMol.java path (which needed MarvinBeans on the classpath) with
the cxcalc CLI shipped in Marvin Desktop Suite. One cxcalc process per source,
batch mode.

  input   Biochemistry/Structures/<source>/inchi.tsv   (external_id, inchi, ...)
  output  Biochemistry/Structures/<source>/pkas/marvin_<ver>.tsv
          columns: external_id, kind, value, tool, tool_version
          value:   ";"-joined <fragment>:<atom>:<pKa>, the three-field
                   atom-indexed form these per-source files keep.

cxcalc invocation:

    cxcalc --ignore-error -i x pka --na 20 --nb 20 <file>

--na/--nb 20 because the widest existing 23.4 row carries 19 tokens.

WHY -i: under --ignore-error cxcalc *silently drops* molecules it cannot parse
AND renumbers the surviving rows, so the default `id` column is a position in
the output, not the input -- one unparseable InChI would shift every id after
it onto the wrong compound. Passing -i makes cxcalc emit `idError[N]` where N
is the true 1-based input line, which is what this script anchors on. The
indices are asserted monotonic and in range before anything is written.

ORDERING is cxcalc's own (by significance: apKa1 is the most significant acidic
value), deliberately not re-sorted by atom index. These are microscopic values
-- see pka_encoding.py: they are per-site predictions on one protonation state,
not a ladder -- so significance order is the only ordering that carries
information. Consumers copy the whole string.

ATOM INDICES are in Marvin's atom space, not ours. Quoting pka_encoding.py:
"Marvin reorders atoms on import, so its indices never described our
structures." They are kept here because this file *is* the provenance record;
nothing in the energy path reads the atom slot.

This script does NOT write compound records. Run Update_Compound_pKas.py for
that -- and note that Compounds.loadPerSourcePkas globs every TSV in pkas/ and
lets the LAST-SORTED file win, so a new marvin_26.1.tsv supersedes marvin_23.4
the next time the cascade runs.
"""

SOURCES = ["ChEBI", "KEGG", "MetaCyc", "Rhea"]

if __name__ == "__main__":
    # Validate arguments BEFORE importing anything or touching the database.
    # These scripts mutate the database, and without this an unknown flag or a
    # mistyped mode was silently ignored and the script ran with its defaults:
    # asking Estimate_Reaction_Reversibility.py for --help rewrote 122 files.
    # Placed above the imports so --help works even where a dependency is
    # missing from the path.
    import argparse as _argparse
    _parser = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter)
    _parser.add_argument(
        "sources", nargs="*", choices=SOURCES, default=SOURCES,
        help="sources to regenerate (default: all four)")
    _ARGS = _parser.parse_args()


import os
import re
import subprocess
import sys
import tempfile
from csv import DictReader

TOOL = "Marvin"
NUM_ACIDIC = 20
NUM_BASIC = 20

# ChEBI pKa rows are keyed CHEBI_<id> while inchi.tsv carries the bare id;
# every other source already agrees. loadPerSourcePkas strips this back off.
ID_PREFIX = {"ChEBI": "CHEBI_"}

STRUCT_ROOT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "Biochemistry", "Structures"
)

ID_RE = re.compile(r"idError\[(\d+)\]")


def marvin_version(cxcalc="cxcalc"):
    """Major.minor of the installed suite, for the filename and the tool_version column."""
    out = subprocess.run(
        [cxcalc, "--help"], capture_output=True, text=True
    ).stdout
    m = re.search(r"(\d+\.\d+)", out)
    return m.group(1) if m else "unknown"


def read_inchis(source):
    """Return [(external_id, inchi)] in file order, skipping rows with no structure."""
    path = os.path.join(STRUCT_ROOT, source, "inchi.tsv")
    rows = []
    with open(path) as fh:
        for line in DictReader(fh, dialect="excel-tab"):
            ext_id, inchi = line.get("external_id"), line.get("inchi")
            if ext_id and inchi:
                rows.append((ext_id, inchi))
    return rows


def run_cxcalc(inchis, cxcalc="cxcalc"):
    """Feed one InChI per line; return the parsed TSV rows as dicts."""
    with tempfile.NamedTemporaryFile("w", suffix=".inchi", delete=False) as fh:
        fh.write("\n".join(inchis) + "\n")
        tmp = fh.name
    try:
        proc = subprocess.run(
            [cxcalc, "--ignore-error", "-i", "x",
             "pka", "--na", str(NUM_ACIDIC), "--nb", str(NUM_BASIC), tmp],
            capture_output=True, text=True,
        )
    finally:
        os.unlink(tmp)
    if not proc.stdout.strip():
        raise RuntimeError(f"empty cxcalc output; stderr:\n{proc.stderr[:2000]}")
    return list(DictReader(proc.stdout.splitlines(), dialect="excel-tab"))


def split_values(row):
    """(acidic, basic, atoms) for one cxcalc row, or None if it carries no usable pKa.

    The atoms column lists one index per REPORTED value, acidic first then
    basic, so it is zipped against acidic+basic in that order.
    """
    if any("FAILED" in (v or "") for v in row.values()):
        return None
    acidic = [row[f"apKa{i}"] for i in range(1, NUM_ACIDIC + 1) if (row.get(f"apKa{i}") or "").strip()]
    basic = [row[f"bpKa{i}"] for i in range(1, NUM_BASIC + 1) if (row.get(f"bpKa{i}") or "").strip()]
    if not acidic and not basic:
        return None
    atoms = [a for a in (row.get("atoms") or "").split(",") if a.strip()]
    if len(atoms) != len(acidic) + len(basic):
        return None                      # misaligned; refuse to guess
    return acidic, basic, atoms


def encode(values, atoms, fragment=1):
    return ";".join(f"{fragment}:{a}:{v}" for v, a in zip(values, atoms))


def process(source, version, cxcalc="cxcalc"):
    rows = read_inchis(source)
    out_rows = []

    results = run_cxcalc([i for _, i in rows], cxcalc=cxcalc)

    # Anchor every result on the true input line, never on output position.
    indices = []
    for r in results:
        m = ID_RE.match((r.get("x") or r.get("id") or "").strip())
        if not m:
            raise RuntimeError(f"{source}: cannot read input index from {r!r}")
        indices.append(int(m.group(1)))
    if indices != sorted(indices) or len(set(indices)) != len(indices):
        raise RuntimeError(f"{source}: cxcalc input indices are not strictly increasing")
    if indices and (min(indices) < 1 or max(indices) > len(rows)):
        raise RuntimeError(f"{source}: cxcalc input index out of range 1..{len(rows)}")

    dropped = len(rows) - len(results)
    prefix = ID_PREFIX.get(source, "")

    for idx, r in zip(indices, results):
        ext_id = prefix + rows[idx - 1][0]
        parsed = split_values(r)
        if parsed is None:
            continue
        acidic, basic, atoms = parsed
        if acidic:
            out_rows.append((ext_id, "pKa", encode(acidic, atoms[: len(acidic)]), TOOL, version))
        if basic:
            out_rows.append((ext_id, "pKb", encode(basic, atoms[len(acidic):]), TOOL, version))

    out_dir = os.path.join(STRUCT_ROOT, source, "pkas")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"marvin_{version}.tsv")
    with open(out_path, "w") as fh:
        fh.write("external_id\tkind\tvalue\ttool\ttool_version\n")
        for row in out_rows:
            fh.write("\t".join(row) + "\n")

    ids = len({r[0] for r in out_rows})
    print(f"{source:<8} in={len(rows):<6} parsed={len(results):<6} unparseable={dropped:<4} "
          f"compounds_with_pka={ids:<6} rows={len(out_rows):<6} -> {out_path}")
    return out_path


def main():
    cxcalc = os.environ.get("CXCALC_BIN", "cxcalc")
    version = marvin_version(cxcalc)
    print(f"Marvin {version} via {cxcalc}\n")
    for source in _ARGS.sources:
        process(source, version, cxcalc=cxcalc)


if __name__ == "__main__":
    main()
