#!/usr/bin/env python
"""Replace Marvin 26.1's pKb on the site classes it gets wrong with 23.4's.

  input   Biochemistry/Structures/<source>/pkas/marvin_<ver>.tsv   (rewritten)
          Biochemistry/Structures/<source>/pkas/marvin_23.4.tsv    (reference)
          Biochemistry/Structures/<source>/{inchi,smiles}.tsv      (atom space)
  report  Biochemistry/Structures/_reports/site_pka_override_<source>.tsv

Site classes and why: see Site_pKa_Overrides.py (adenine N1: 26.1 7.0-7.3 vs
Alberty 4.0-4.7; histidine-type imidazole: 26.1 7.20-7.30 vs IUPAC 6.04-6.05).
The ladders feed the thermodynamic pH transforms, so they must agree with the
pH-7 structures, which Run_Marvin_Protonations.py corrects for the same classes.

ATOM SPACE. The ladder's atom slot is Marvin's own import order of the
structure Run_Marvin_pKas.py fed it: for a compound with an InChI, the
canonical SMILES Run_Marvin_Protonations.pka_input() derives from it (cxcalc
CLI, fragment always 1; the raw InChI if RDKit cannot read it), else the
source SMILES (PkaPlugin, one fragment number per convertToFrags() piece).
This script imports the same string the same way and locates each class site
in that order. The 23.4 ladders were produced in a different atom space that
cannot be rebuilt, so 23.4's value is found by value instead: the single basic
value in the class window (Site_pKa_Overrides.CLASSES[...]["window_234"]).
Where 23.4 has none, or more than one, the median of the unambiguous 23.4
values over ALL FOUR sources is used -- whichever sources a run rewrites --
so a partial rerun writes the same numbers as a full one.

Only a site token at or above the class threshold is replaced; the basic
tokens are then re-sorted descending, the significance order the ladders use.
Idempotent: a second run finds the sites below threshold and changes nothing,
and the report of an earlier run is kept rather than overwritten with nothing.

Usage:  Override_Site_pKas.py [sources...] [--version 26.1]
"""

import argparse
import csv
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import Site_pKa_Overrides as SPO                               # noqa: E402
from Site_pKa_Overrides import (CLASSES, format_tokens, median_234,  # noqa: E402
                                parse_tokens, value_234)
from rdkit import Chem, RDLogger                               # noqa: E402

RDLogger.DisableLog("rdApp.*")

SOURCES = ["ChEBI", "KEGG", "MetaCyc", "Rhea"]
ID_PREFIX = {"ChEBI": "CHEBI_"}            # as Run_Marvin_pKas.py writes ids
STRUCT_ROOT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "Biochemistry", "Structures")
REPORT_HEADER = "external_id\tclass\tfragment\tatom\twas_pkb\tnow_pkb\tbasis\n"


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
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


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


def _classes_in(struct):
    """Site classes a structure string contains (RDKit, no Marvin)."""
    m = (Chem.MolFromInchi(struct) if struct.startswith("InChI=")
         else Chem.MolFromSmiles(struct, sanitize=False))
    if m is None:
        return set()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return {name for name, _ in SPO.find_sites(m)}


def site_slots(struct, path, MolImporter, MolExporter):
    """{(fragment, atom): class} ladder slots of every class site, as strings."""
    if path == "cxcalc":
        # Run_Marvin_pKas.py feeds cxcalc the canonical-tautomer SMILES, not
        # the InChI; import the identical string so the atom order matches.
        from Run_Marvin_Protonations import pka_input
        struct = pka_input(struct) or struct
    mol = MolImporter.importMol(struct)
    frags = [mol] if path == "cxcalc" else list(mol.convertToFrags())
    out = {}
    for fi, frag in enumerate(frags, start=1):
        for name, site in SPO.find_sites(SPO.rdkit_view(frag, MolExporter)):
            for i in site:
                out[(str(fi), str(i + 1))] = name
    return out


def collect(source, version):
    """Per-source work list: class-bearing ids with their 26.1 rows and 23.4 pKb ladders."""
    pdir = os.path.join(STRUCT_ROOT, source, "pkas")
    new_path = os.path.join(pdir, f"marvin_{version}.tsv")
    old_path = os.path.join(pdir, "marvin_23.4.tsv")
    rows = _read_ladders(new_path) if os.path.exists(new_path) else []
    old = ({r["external_id"]: r["value"] for r in _read_ladders(old_path) if r["kind"] == "pKb"}
           if os.path.exists(old_path) else {})
    structs = _structures(source)
    members = {}
    for r in rows:
        if r["kind"] != "pKb":
            continue
        s = structs.get(r["external_id"])
        if s:
            cl = _classes_in(s[0])
            if cl:
                members[r["external_id"]] = (s, cl)
    return new_path, rows, old, members


def apply(source, version, fallbacks, MolImporter, MolExporter, work):
    new_path, rows, old, members = work
    report = []
    for r in rows:
        eid = r["external_id"]
        if r["kind"] != "pKb" or eid not in members:
            continue
        (struct, path), _ = members[eid]
        slots = site_slots(struct, path, MolImporter, MolExporter)
        toks = parse_tokens(r["value"])
        changed = False
        for k, (f, a, v) in enumerate(toks):
            name = slots.get((f, a))
            if name is None or v < CLASSES[name]["threshold_261"]:
                continue
            v234 = value_234(name, old.get(eid, ""))
            new_val, basis = ((v234, "23.4_same_compound") if v234 is not None
                              else (fallbacks[name], "23.4_median"))
            toks[k] = (f, a, new_val)
            report.append((eid, name, f, a, f"{v:.2f}", f"{new_val:.2f}", basis))
            changed = True
        if changed:
            toks.sort(key=lambda t: -t[2])
            r["value"] = format_tokens(toks)

    if rows:
        with open(new_path, "w") as fh:
            fh.write("external_id\tkind\tvalue\ttool\ttool_version\n")
            for r in rows:
                fh.write("\t".join([r["external_id"], r["kind"], r["value"], r["tool"],
                                    r["tool_version"]]) + "\n")
    rep_dir = os.path.join(STRUCT_ROOT, "_reports")
    os.makedirs(rep_dir, exist_ok=True)
    rep_path = os.path.join(rep_dir, f"site_pka_override_{source}.tsv")
    if report or not os.path.exists(rep_path):
        with open(rep_path, "w") as fh:
            fh.write(REPORT_HEADER)
            for row in report:
                fh.write("\t".join(row) + "\n")
        kept = ""
    else:
        kept = "  (nothing replaced; earlier report kept)"
    by = {n: sum(1 for x in report if x[1] == n) for n in CLASSES}
    print(f"{source:<8} class-bearing ladders={len(members):<5} replaced: "
          + " ".join(f"{n}={c}" for n, c in by.items())
          + f"  same_compound={sum(1 for x in report if x[6] == '23.4_same_compound')}"
          + f" median={sum(1 for x in report if x[6] == '23.4_median')}{kept}")
    return report


def run(sources, version):
    MolImporter, MolExporter = _marvin()
    # Fallback medians pool ALL FOUR sources' unambiguous 23.4 values, so a
    # partial rerun writes the same numbers as a full one.
    work_all = {s: collect(s, version) for s in SOURCES}
    fallbacks = {}
    for name in CLASSES:
        pooled = [w[2][e] for w in work_all.values() for e, (_, cl) in w[3].items()
                  if name in cl and e in w[2]]
        fallbacks[name] = median_234(name, pooled)
    print("23.4 fallback medians: " + ", ".join(f"{n}={v}" for n, v in fallbacks.items()))
    for s in sources:
        apply(s, version, fallbacks, MolImporter, MolExporter, work_all[s])


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
