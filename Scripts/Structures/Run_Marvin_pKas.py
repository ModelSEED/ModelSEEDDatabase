#!/usr/bin/env python
"""Regenerate the per-source Marvin pKa bundles with ChemAxon cxcalc.

Replaces the pKaMol.java path (which needed MarvinBeans on the classpath) with
the cxcalc CLI shipped in Marvin Desktop Suite. One cxcalc process per input
representation per source, batch mode.

  input   Biochemistry/Structures/<source>/inchi.tsv    (preferred)
          Biochemistry/Structures/<source>/smiles.tsv   (gap-fill only)
  output  Biochemistry/Structures/<source>/pkas/marvin_<ver>.tsv
          columns: external_id, kind, value, tool, tool_version
          value:   ";"-joined <fragment>:<atom>:<pKa>, the three-field
                   atom-indexed form these per-source files keep.

INPUT SELECTION. `smiles.tsv` is a strict superset of `inchi.tsv`: every InChI
id also has a SMILES, and 8,834 compounds across the four sources have a SMILES
and no InChI. Those 8,834 are why the 23.4 bundles cover more ids than an
InChI-only run.

The cxcalc CLI cannot process any of them, and the reason they have no InChI is
the same reason cxcalc refuses them:

  8,728  carry `*` attachment points (structural repeating units). Marvin reads
         these as QUERY molecules and cxcalc declines outright:
         "pka: Calculation result is not defined for query molecules".
    105  are organometallics (Mg-porphyrins, chlorophylls). They parse, but
         cxcalc returns no pKa and an empty `atoms` column.
      1  is a dative-bond SMILES ("...[Mg]35<-N2=...") the SMILES parser
         rejects at the '<' character.

The Java API underneath the CLI has no such restriction. So those compounds go
through chemaxon.calculations.PkaPlugin directly, via JPype (a JRE is enough --
no JDK, nothing to compile). This is the same route pKaMol.java took before the
API was renamed in 26.1 and its committed .class stopped running.

Both paths are ONE engine, verified rather than assumed: over 300 ChEBI
compounds the plugin reproduces cxcalc to a median |delta| of 0.0000 and a max
of 0.000, 100% within 0.01, with 1 site-count difference in 567 sets and
identical atom ordering in 562. A bundle may therefore mix them. Which rows
came from the plugin is recoverable from the data: they are exactly the ids in
the bundle that are absent from `inchi.tsv`.

The plugin does NOT invent sites on the wildcards: across the Rhea polymers,
0 of 154 predicted sites sat on an atom bonded to a `*`, and site counts match
23.4 in 40 of 40 comparable sets.

InChI and SMILES are also NOT interchangeable where both exist. Measured on 600
ChEBI compounds, they agree within 0.05 for 80.7% of sites when the SMILES is
neutral, but only 53.1% when it carries a charge, with a 27% site-count
mismatch: a charged SMILES is an already-deprotonated species, so Marvin is
answering a different question than it is for the neutral InChI parent. InChI
therefore wins wherever it exists, and no InChI-derived value is displaced.

REQUIRES. `cxcalc` on PATH (Marvin Desktop Suite). For the plugin path only,
`pip install jpype1`, plus the suite's jars -- found next to the cxcalc binary,
or point MARVIN_LIB at <marvinsuite>/lib. `--structures inchi` skips the plugin
path entirely and needs neither.

cxcalc invocation:

    cxcalc --ignore-error -i x pka --na 20 --nb 20 <file>

--na/--nb 20 because the widest existing 23.4 row carries 19 tokens.

WHY -i: under --ignore-error cxcalc *silently drops* molecules it cannot parse
AND renumbers the surviving rows, so the default `id` column is a position in
the output, not the input -- one unparseable structure would shift every id
after it onto the wrong compound. Passing -i makes cxcalc emit `idError[N]`
where N is the true 1-based input line, which is what this script anchors on.
The indices are asserted monotonic and in range before anything is written.

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
    # No choices= here: with nargs="*" argparse validates the DEFAULT against
    # choices as a single value, so any default fails. Validated by hand below.
    _parser.add_argument(
        "sources", nargs="*", default=[],
        help=f"sources to regenerate, from {SOURCES} (default: all four)")
    _parser.add_argument(
        "--structures", choices=["inchi", "smiles", "both"], default="both",
        help="input representation: 'inchi' only, 'smiles' only, or 'both' "
             "(default) = InChI where it exists, SMILES to fill the gap")
    _ARGS = _parser.parse_args()
    _bad = [s for s in _ARGS.sources if s not in SOURCES]
    if _bad:
        _parser.error(f"invalid source(s) {_bad}; choose from {SOURCES}")


import glob
import os
import re
import shutil
import subprocess
import tempfile
from csv import DictReader

TOOL = "Marvin"
NUM_ACIDIC = 20
NUM_BASIC = 20

# ChEBI pKa rows are keyed CHEBI_<id> while the structure files carry the bare
# id; every other source already agrees. loadPerSourcePkas strips this back off.
ID_PREFIX = {"ChEBI": "CHEBI_"}

STRUCT_ROOT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "Biochemistry", "Structures"
)

ID_RE = re.compile(r"idError\[(\d+)\]")


def marvin_version(cxcalc="cxcalc"):
    """Major.minor of the installed suite, for the filename and tool_version column."""
    out = subprocess.run([cxcalc, "--help"], capture_output=True, text=True).stdout
    m = re.search(r"(\d+\.\d+)", out)
    return m.group(1) if m else "unknown"


def read_structures(source, column, filename):
    """Return [(external_id, structure)] in file order, skipping empty rows."""
    path = os.path.join(STRUCT_ROOT, source, filename)
    rows = []
    with open(path) as fh:
        for line in DictReader(fh, dialect="excel-tab"):
            ext_id, struct = line.get("external_id"), line.get(column)
            if ext_id and struct:
                rows.append((ext_id, struct))
    return rows


def run_cxcalc(structures, suffix, cxcalc="cxcalc"):
    """Feed one structure per line; return (input_index, row) pairs.

    The index is the TRUE 1-based input line, read back from cxcalc's
    idError[N] id column -- never the output row position, which renumbers
    whenever --ignore-error drops an unparseable structure.
    """
    if not structures:
        return []
    with tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False) as fh:
        fh.write("\n".join(structures) + "\n")
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
        # Every structure in the batch failed to parse. That is a legitimate
        # outcome for a gap-fill batch (see INPUT SELECTION: the SMILES-only
        # compounds are query molecules and organometallics), so report it
        # rather than aborting the source. An empty stdout with an empty
        # stderr is different -- that means cxcalc itself failed.
        if proc.stderr.strip():
            print(f"  ! cxcalc returned nothing for all {len(structures)} structures "
                  f"in this batch; first error: {proc.stderr.splitlines()[0][:120]}")
            return []
        raise RuntimeError("cxcalc produced no output and no error")

    out = []
    for r in DictReader(proc.stdout.splitlines(), dialect="excel-tab"):
        m = ID_RE.match((r.get("x") or r.get("id") or "").strip())
        if not m:
            raise RuntimeError(f"cannot read input index from {r!r}")
        out.append((int(m.group(1)), r))

    idx = [i for i, _ in out]
    if idx != sorted(idx) or len(set(idx)) != len(idx):
        raise RuntimeError("cxcalc input indices are not strictly increasing")
    if idx and (min(idx) < 1 or max(idx) > len(structures)):
        raise RuntimeError(f"cxcalc input index out of range 1..{len(structures)}")
    return out


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


def _marvin_lib(cxcalc="cxcalc"):
    """The jar directory of the installed suite, for the plugin classpath."""
    lib = os.environ.get("MARVIN_LIB")
    if lib:
        return lib
    exe = shutil.which(cxcalc)
    if exe:
        cand = os.path.join(os.path.dirname(os.path.realpath(exe)), "..", "lib")
        if os.path.isdir(cand):
            return os.path.normpath(cand)
    raise RuntimeError(
        "cannot locate the Marvin jars; set MARVIN_LIB to <marvinsuite>/lib")


def run_plugin(structures, cxcalc="cxcalc"):
    """pKa via the Java PkaPlugin, for structures the cxcalc CLI refuses.

    Returns [(input_index, {fragment: (acidic, basic)})] where each of acidic
    and basic is [(atom, value)] -- acidic ascending, basic descending, which
    is the significance order cxcalc itself emits.

    jpype is imported lazily so --help and the InChI-only path never start a
    JVM. A JRE is sufficient; nothing here is compiled.
    """
    if not structures:
        return []
    import jpype                                   # noqa: local by design
    import jpype.imports                           # noqa: registers importer
    if not jpype.isJVMStarted():
        jpype.startJVM(classpath=glob.glob(os.path.join(_marvin_lib(cxcalc), "*.jar")))
    from chemaxon.calculations import PkaPlugin
    from chemaxon.formats import MolImporter

    out = []
    for i, struct in enumerate(structures, start=1):
        try:
            mol = MolImporter.importMol(struct)
            frags = {}
            for fi, frag in enumerate(mol.convertToFrags(), start=1):
                plugin = PkaPlugin()
                plugin.setpH(7.0)
                plugin.setMolecule(frag)
                plugin.run()
                acidic, basic = [], []
                for atom in range(frag.getAtomCount()):
                    a = plugin.getpKaValues(atom, PkaPlugin.ACIDIC)
                    b = plugin.getpKaValues(atom, PkaPlugin.BASIC)
                    if a is not None:
                        acidic.append((atom + 1, round(float(a[0]), 2)))
                    if b is not None:
                        basic.append((atom + 1, round(float(b[0]), 2)))
                acidic.sort(key=lambda t: t[1])
                basic.sort(key=lambda t: -t[1])
                if acidic or basic:
                    frags[fi] = (acidic, basic)
            if frags:
                out.append((i, frags))
        except Exception:
            continue                               # unreadable structure
    return out


def emit_plugin(rows, results, prefix, version):
    """Plugin results -> output tuples, tokens spanning fragments as 23.4 did."""
    out = []
    for idx, frags in results:
        ext_id = prefix + rows[idx - 1][0]
        for kind, which in (("pKa", 0), ("pKb", 1)):
            toks = [f"{fi}:{atom}:{val}"
                    for fi in sorted(frags)
                    for atom, val in frags[fi][which]]
            if toks:
                out.append((ext_id, kind, ";".join(toks), TOOL, version))
    return out


def emit(rows, results, prefix, version):
    """Turn (index, cxcalc row) pairs into output tuples."""
    out = []
    for idx, r in results:
        parsed = split_values(r)
        if parsed is None:
            continue
        acidic, basic, atoms = parsed
        ext_id = prefix + rows[idx - 1][0]
        if acidic:
            out.append((ext_id, "pKa", encode(acidic, atoms[: len(acidic)]), TOOL, version))
        if basic:
            out.append((ext_id, "pKb", encode(basic, atoms[len(acidic):]), TOOL, version))
    return out


def process(source, version, structures="both", cxcalc="cxcalc"):
    prefix = ID_PREFIX.get(source, "")
    out_rows = []
    stats = {}

    inchi_rows, smiles_rows = [], []
    if structures in ("inchi", "both"):
        inchi_rows = read_structures(source, "inchi", "inchi.tsv")
    if structures in ("smiles", "both"):
        smiles_rows = read_structures(source, "smiles", "smiles.tsv")
        if structures == "both":
            # InChI wins wherever it exists; SMILES fills the gap only.
            have = {e for e, _ in inchi_rows}
            smiles_rows = [(e, s) for e, s in smiles_rows if e not in have]

    # InChI through the cxcalc CLI; the SMILES-only remainder through the Java
    # PkaPlugin, which accepts the query molecules and organometallics the CLI
    # refuses. Same engine either way -- see INPUT SELECTION.
    if inchi_rows:
        results = run_cxcalc([s for _, s in inchi_rows], ".inchi", cxcalc=cxcalc)
        got = emit(inchi_rows, results, prefix, version)
        out_rows += got
        stats["inchi (cxcalc)"] = (
            len(inchi_rows), len(inchi_rows) - len(results), len({r[0] for r in got}))

    if smiles_rows:
        results = run_plugin([s for _, s in smiles_rows], cxcalc=cxcalc)
        got = emit_plugin(smiles_rows, results, prefix, version)
        out_rows += got
        # "no_result" is not all failure: of the 191 across all four sources,
        # 1 SMILES is unparseable and 1 crashes PkaPlugin.run(); the other 189
        # run cleanly and simply have no ionizable site.
        stats["smiles (plugin)"] = (
            len(smiles_rows), len(smiles_rows) - len(results), len({r[0] for r in got}))

    out_dir = os.path.join(STRUCT_ROOT, source, "pkas")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"marvin_{version}.tsv")
    with open(out_path, "w") as fh:
        fh.write("external_id\tkind\tvalue\ttool\ttool_version\n")
        for row in out_rows:
            fh.write("\t".join(row) + "\n")

    parts = "  ".join(
        f"{k}: in={v[0]} no_result={v[1]} with_pka={v[2]}" for k, v in stats.items())
    print(f"{source:<8} {parts}")
    print(f"{'':<8} total compounds={len({r[0] for r in out_rows}):<6} "
          f"rows={len(out_rows):<6} -> {os.path.relpath(out_path, STRUCT_ROOT)}")
    return out_path


def main():
    cxcalc = os.environ.get("CXCALC_BIN", "cxcalc")
    version = marvin_version(cxcalc)
    print(f"Marvin {version} via {cxcalc}  (structures={_ARGS.structures})\n")
    for source in (_ARGS.sources or SOURCES):
        process(source, version, structures=_ARGS.structures, cxcalc=cxcalc)


if __name__ == "__main__":
    main()
