#!/usr/bin/env python
"""Replace Marvin 26.1's adenine N1 pKb in the pKa ladders with 23.4's.

  input   Biochemistry/Structures/<source>/pkas/marvin_<ver>.tsv   (rewritten)
          Biochemistry/Structures/<source>/pkas/marvin_23.4.tsv    (reference)
          Biochemistry/Structures/<source>/{inchi,smiles}.tsv      (atom space)
  report  Biochemistry/Structures/_reports/adenine_n1_pka_override_<source>.tsv

Why: see Adenine_N1.py. 26.1 puts adenine N1 at 7.0-7.3, 23.4 at 2.2-2.5
(measured ~3.5-4.2). The ladders feed the thermodynamic pH transforms, so
leaving 26.1's value there would make them disagree with the protonation
bundles, which Run_Marvin_Protonations.py already corrects.

ATOM SPACE. The ladder's atom slot is Marvin's own import order of the
structure Run_Marvin_pKas.py fed it: for a compound with an InChI, the
canonical-tautomer SMILES Run_Marvin_Protonations.pka_input() derives from it
(cxcalc CLI, fragment always 1; the raw InChI if RDKit cannot read it), else
the source SMILES (PkaPlugin, one fragment number per convertToFrags() piece). This script imports the same string the same way and
locates N1 in that order, so the token it rewrites is the one on N1. The
23.4 ladders were produced in a different atom space, so their N1 is found by
value instead: the single basic value in Adenine_N1.N1_234_RANGE. Where 23.4
has none, or more than one, the median of the unambiguous 23.4 values is used,
and the report says which. (The 23.4 atom space cannot be rebuilt: no
import order of the current InChI or SMILES puts N1 on 23.4's atom slot.)
The common ambiguous case is the acyl-CoAs, whose 23.4 ladders carry 1.74
and 2.87. Which one is N1 does not matter downstream: any N1 value from 1.7
to 2.9 leaves the site >99.9% unprotonated at pH 7, and changes the pH
transform of the formation energy by <0.01 kJ/mol.

Only an N1 token at or above Adenine_N1.N1_26_THRESHOLD is replaced; the
basic tokens are then re-sorted descending, the significance order the
ladders already use. Idempotent: a second run finds N1 below the threshold
and changes nothing.

Usage:  Override_Adenine_N1_pKas.py [sources...] [--version 26.1]
"""

import argparse
import csv
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import Adenine_N1                                              # noqa: E402
from Adenine_N1 import (format_tokens, median_n1_234, n1_value_234,  # noqa: E402
                        parse_tokens, N1_26_THRESHOLD)
from rdkit import Chem, RDLogger                               # noqa: E402

RDLogger.DisableLog("rdApp.*")

SOURCES = ["ChEBI", "KEGG", "MetaCyc", "Rhea"]
ID_PREFIX = {"ChEBI": "CHEBI_"}            # as Run_Marvin_pKas.py writes ids
STRUCT_ROOT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "Biochemistry", "Structures")


def _marvin():
    import jpype
    if not jpype.isJVMStarted():
        lib = os.environ.get("MARVIN_LIB")
        if not lib:
            exe = shutil.which(os.environ.get("CXCALC_BIN", "cxcalc"))
            if not exe:
                raise SystemExit("cannot locate the Marvin jars; set MARVIN_LIB")
            lib = os.path.join(os.path.dirname(os.path.realpath(exe)), "..", "lib")
        jpype.startJVM(classpath=glob.glob(os.path.join(lib, "*.jar")))
    return (jpype.JClass("chemaxon.formats.MolImporter"),
            jpype.JClass("chemaxon.formats.MolExporter"))


def _read_ladders(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    return rows


def _structures(source):
    """{ext_id (with prefix): (structure, path)} in Run_Marvin_pKas.py's choice."""
    prefix = ID_PREFIX.get(source, "")
    out = {}
    for fname, col, path in (("smiles.tsv", "smiles", "plugin"), ("inchi.tsv", "inchi", "cxcalc")):
        p = os.path.join(STRUCT_ROOT, source, fname)
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get(col):
                    out[prefix + r["external_id"]] = (r[col], path)   # InChI overwrites
    return out


def _has_adenine(struct):
    m = (Chem.MolFromInchi(struct) if struct.startswith("InChI=")
         else Chem.MolFromSmiles(struct, sanitize=False))
    if m is None:
        return False
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m.HasSubstructMatch(Adenine_N1.ADENINE)


def n1_tokens(struct, path, MolImporter, MolExporter):
    """[(fragment, atom)] ladder slots of every adenine N1, as strings."""
    if path == "cxcalc":
        # Run_Marvin_pKas.py feeds cxcalc the canonical-tautomer SMILES, not
        # the InChI; import the identical string so the atom order matches.
        from Run_Marvin_Protonations import pka_input
        struct = pka_input(struct) or struct
    mol = MolImporter.importMol(struct)
    if path == "cxcalc":
        return [("1", str(i + 1)) for i in Adenine_N1.n1_atoms(mol, MolExporter)]
    out = []
    for fi, frag in enumerate(mol.convertToFrags(), start=1):
        out += [(str(fi), str(i + 1)) for i in Adenine_N1.n1_atoms(frag, MolExporter)]
    return out


def collect(source, version):
    """Per-source work list: adenine ids with their 26.1 and 23.4 pKb ladders."""
    pdir = os.path.join(STRUCT_ROOT, source, "pkas")
    new_path = os.path.join(pdir, f"marvin_{version}.tsv")
    old_path = os.path.join(pdir, "marvin_23.4.tsv")
    rows = _read_ladders(new_path)
    old = {r["external_id"]: r["value"] for r in _read_ladders(old_path)
           if r["kind"] == "pKb"} if os.path.exists(old_path) else {}
    structs = _structures(source)
    adenine = {}
    for r in rows:
        if r["kind"] != "pKb":
            continue
        s = structs.get(r["external_id"])
        if s and _has_adenine(s[0]):
            adenine[r["external_id"]] = s
    return new_path, rows, old, adenine


def apply(source, version, fallback, MolImporter, MolExporter, work):
    new_path, rows, old, adenine = work
    report = []
    for r in rows:
        eid = r["external_id"]
        if r["kind"] != "pKb" or eid not in adenine:
            continue
        struct, path = adenine[eid]
        slots = set(n1_tokens(struct, path, MolImporter, MolExporter))
        toks = parse_tokens(r["value"])
        v234 = n1_value_234(old.get(eid, ""))
        new_val, basis = (v234, "23.4_same_compound") if v234 is not None else (fallback, "23.4_median")
        changed = False
        for k, (f, a, v) in enumerate(toks):
            if (f, a) in slots and v >= N1_26_THRESHOLD:
                toks[k] = (f, a, new_val)
                report.append((eid, f, a, f"{v:.2f}", f"{new_val:.2f}", basis))
                changed = True
        if not any((f, a) in slots for f, a, _ in toks):
            report.append((eid, "", "", "", "", "no_n1_token"))
        if changed:
            toks.sort(key=lambda t: -t[2])
            r["value"] = format_tokens(toks)

    with open(new_path, "w") as fh:
        fh.write("external_id\tkind\tvalue\ttool\ttool_version\n")
        for r in rows:
            fh.write("\t".join([r["external_id"], r["kind"], r["value"], r["tool"],
                                r["tool_version"]]) + "\n")
    rep_dir = os.path.join(STRUCT_ROOT, "_reports")
    os.makedirs(rep_dir, exist_ok=True)
    with open(os.path.join(rep_dir, f"adenine_n1_pka_override_{source}.tsv"), "w") as fh:
        fh.write("external_id\tfragment\tatom\twas_pkb\tnow_pkb\tbasis\n")
        for row in report:
            fh.write("\t".join(row) + "\n")
    n = sum(1 for x in report if x[5] != "no_n1_token")
    print(f"{source:<8} adenine ladders={len(adenine):<5} N1 replaced={n:<5} "
          f"same_compound={sum(1 for x in report if x[5] == '23.4_same_compound'):<5} "
          f"median={sum(1 for x in report if x[5] == '23.4_median'):<4} "
          f"no_n1_token={sum(1 for x in report if x[5] == 'no_n1_token')}")
    return report


def run(sources, version):
    MolImporter, MolExporter = _marvin()
    work = {s: collect(s, version) for s in sources}
    # The fallback is pooled over every source's unambiguous adenine ladders.
    pooled = {}
    for s, (_, _, old, adenine) in work.items():
        pooled.update({f"{s}:{e}": old[e] for e in adenine if e in old})
    fallback = median_n1_234(pooled)
    print(f"23.4 adenine N1 median (fallback) = {fallback}")
    for s in sources:
        apply(s, version, fallback, MolImporter, MolExporter, work[s])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sources", nargs="*", default=[])
    ap.add_argument("--version", default="26.1")
    a = ap.parse_args()
    bad = [s for s in a.sources if s not in SOURCES]
    if bad:
        ap.error(f"invalid source(s) {bad}; choose from {SOURCES}")
    run(a.sources or SOURCES, a.version)


if __name__ == "__main__":
    main()
