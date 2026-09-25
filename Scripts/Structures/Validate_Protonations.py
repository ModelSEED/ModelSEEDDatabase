#!/usr/bin/env python
"""Check every protonation row against the chemistry a protonation can do.

WHY. A protonation moves protons and nothing else. Between a source structure
and its pH-adjusted form, therefore, the change in hydrogen count must equal
the change in charge, and no other element may change:

    dH == dcharge        and        heavy atoms conserved

Any row that breaks this did not undergo a protonation. It underwent something
else -- and in this repository "something else" has one recurring cause: InChI
disconnects metal-ligand bonds by design, so Marvin is handed the ligands as
free ions and protonates them on their own. Fe4S4 comes back as H8Fe4S4 (eight
hydrogens gained, charge unchanged), triphenyltin chloride as C18H18ClSn, and
heptamolybdate as H48Mo7O24. Bare S, Se, P and O atoms come back as H2S, H2Se,
PH3 and H2O in BOTH representations. None of these is a protonation state; all
of them ship as a compound's formula if nothing checks.

This script is that check. It is representation-agnostic and needs no list of
metals, which is the point: whatever future engine change or source quirk
produces an impossible row, the invariant catches it.

THREE CHECKS, in decreasing order of what they mean for the release:

  row       a protonation row against its own representation's source row
            (InChI row vs inchi.tsv, SMILE row vs smiles.tsv). A failure is a
            row that must not feed a compound's formula.

  bundle    the InChI row against the SMILE row of the same compound in the
            same bundle. They describe the same molecule at the same pH, so
            they must agree by the same rule. A failure here with both rows
            individually passing means the two representations disagree about
            what the molecule is -- usually the source does too (next check).

  layers    what an inchi.tsv row STORES against what its InChI string
            DECLARES in its own formula, /p and /q layers. These must agree
            by definition; when they do not, the stored column came from a
            parser that failed the string. Found because ChEBI stored
            chlorate, InChI=1S/ClO3/c2-1(3)4/q-1, with charge 0 (RDKit
            rejects the hypervalent chlorine, the OpenBabel fallback drops
            the charge with a warning), and RDKit stored the Mg porphyrins
            one proton and one charge unit too high. Print_Structure_
            Formula_Charge.parse_structure now derives InChI rows from the
            layers, so this check is the regression guard for that.

  key       the bundle's InChIKey row against the key of its InChI row. The
            key is a hash of the string, so they must agree; when they do
            not, one of the two was rewritten without the other. Added after
            the first post-hoc repair replaced InChI rows and left their
            keys hashed from the discarded protonated string -- elemental
            sulfur came out keyed as hydrosulfide.

  source    inchi.tsv against smiles.tsv for the same compound. The two
            representations may differ in protonation state (the InChI is the
            neutral parent, the SMILES the charged species) but not in heavy
            atoms, and their H and charge must differ together. MetaCyc CPD-7
            fails this: Fe4S4 in InChI, H4Fe4S4 in SMILES with explicit [SH]
            on bridging sulfides. That is a curation question, not a Marvin
            one, and it is why "prefer the SMILE row" cannot fix everything.

For every row failure the report says whether the OTHER representation's row
passes, because that is the cheapest fix available: prefer the row that obeys
chemistry. Compounds are resolved through the alias file so a failure can be
read as a ModelSEED id and counted against reactions.

USAGE
    Validate_Protonations.py                         # 26.1, all sources, summary
    Validate_Protonations.py --bundle marvin_23.4_ph7
    Validate_Protonations.py --tsv violations.tsv    # one row per failure
    Validate_Protonations.py --fail-on-violation     # exit 1 on any ROW failure
                                                     # not listed in the allowlist

ALLOWLIST. Biochemistry/Curation/exclusions/protonation_invariant_excluded.tsv
(columns: source, external_id, check, reason, date, curator). A row listed
there is reported but does not fail the run. It exists so a known, curated
case can be tolerated deliberately rather than by everyone learning to ignore
the check.
"""
import argparse
import csv
import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
STRUCT = os.path.join(REPO, "Biochemistry", "Structures")   # overridable with --root
ALIASES = os.path.join(REPO, "Biochemistry", "Aliases", "Unique_ModelSEED_Compound_Aliases.txt")
ALLOWLIST = os.path.join(REPO, "Biochemistry", "Curation", "exclusions",
                         "protonation_invariant_excluded.tsv")
SOURCES = ("MetaCyc", "KEGG", "ChEBI", "Rhea")
ORIGINAL = {"InChI": ("inchi.tsv", "inchi"), "SMILE": ("smiles.tsv", "smiles")}


def parse_formula(f):
    """Element -> count. R groups count as an element so they are conserved too."""
    if not f or f in ("null", "None"):
        return None
    d = defaultdict(int)
    for el, n in re.findall(r"([A-Z][a-z]?|R)(\d*)", f):
        d[el] += int(n) if n else 1
    return dict(d)


def as_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


# Severity of a kind. STRICT kinds are chemistry violations and fail a gated
# run. INFO kinds are reported but never fail it:
#   unchecked   the source row carries no formula/charge, so there is nothing
#               to compare against (108 SMILE rows in 26.1, all Mg porphyrins
#               and the like whose source formula was never derived).
#   wildcard    the structure carries R groups. The source formula column
#               omits R for many of these while the bundle adds it (that is
#               the run script's own R invariant), and the hydrogen count on
#               an attachment atom is convention-dependent. Neither is a
#               protonation error, so the row is reported under this label
#               rather than under a chemistry kind.
INFO_KINDS = ("unchecked", "wildcard")


def compare(f_from, q_from, f_to, q_to):
    """Return (kind, dH, dq). kind is 'ok', an INFO kind, or the chemistry break."""
    P0, P1, Q0, Q1 = parse_formula(f_from), parse_formula(f_to), as_int(q_from), as_int(q_to)
    if P0 is None or P1 is None or Q0 is None or Q1 is None:
        return ("unchecked", None, None)
    dH = P1.get("H", 0) - P0.get("H", 0)
    dq = Q1 - Q0
    heavy = sorted(e for e in set(P0) | set(P1) if e not in ("H", "R") and P0.get(e, 0) != P1.get(e, 0))
    if heavy:
        return ("heavy atoms changed: " + ",".join(heavy), dH, dq)
    if "R" in P0 or "R" in P1:
        return ("ok" if dH == dq else "wildcard", dH, dq)
    if dH == dq:
        return ("ok", dH, dq)
    if dq == 0 and dH > 0:
        return ("H gained, no charge change", dH, dq)
    if dq == 0 and dH < 0:
        return ("H lost, no charge change", dH, dq)
    if dH == 0:
        return ("charge changed, no H change", dH, dq)
    return ("dH != dcharge", dH, dq)


# One implementation of "what does this InChI declare", shared with the
# script that derives the stored columns, so the check and the derivation
# cannot drift apart.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Print_Structure_Formula_Charge import inchi_layers   # noqa: E402


def load_source(src):
    rows = {}
    for typ, (fn, col) in ORIGINAL.items():
        path = os.path.join(STRUCT, src, fn)
        if not os.path.isfile(path):
            continue
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows[(r["external_id"], typ)] = (r["formula"], r["charge"], r.get(col, ""))
    return rows


def load_bundle(src, bundle):
    """(external_id, type) -> (formula, charge, structure); InChIKey rows carry
    the key in the structure slot and empty formula/charge."""
    path = os.path.join(STRUCT, src, "protonations", bundle + ".tsv")
    if not os.path.isfile(path):
        return None
    rows = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["type"] in ORIGINAL and r["formula"]:
                rows[(r["external_id"], r["type"])] = (r["formula"], r["charge"], r["structure"])
            elif r["type"] == "InChIKey" and r["structure"]:
                rows[(r["external_id"], "InChIKey")] = ("", "", r["structure"])
    return rows


def inchikey_of(inchi):
    """Standard InChIKey of an InChI string, or None if it cannot be hashed."""
    try:
        from rdkit import Chem
        return Chem.InchiToInchiKey(inchi) or None
    except Exception:
        return None


def load_aliases():
    """(source, external_id) -> [cpd ids]."""
    m = defaultdict(list)
    if not os.path.isfile(ALIASES):
        return m
    with open(ALIASES) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Source"] in SOURCES:
                m[(r["Source"], r["External ID"])].append(r["ModelSEED ID"])
    return m


def load_allowlist():
    allowed = set()
    if not os.path.isfile(ALLOWLIST):
        return allowed
    with open(ALLOWLIST) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            allowed.add((r.get("source", "").strip(), r.get("external_id", "").strip(),
                         r.get("check", "").strip()))
    return allowed


def reactions_by_compound():
    """cpd -> number of live reactions that use it. Only read when asked for."""
    n = Counter()
    for f in sorted(glob.glob(os.path.join(REPO, "Biochemistry", "reaction_[0-9][0-9].json"))):
        for r in json.load(open(f)):
            if r.get("is_obsolete") == 1:
                continue
            for p in r.get("stoichiometry") or []:
                n[p["compound"]] += 1
    return n


def run(bundle, sources, want_reactions):
    aliases = load_aliases()
    rxn = reactions_by_compound() if want_reactions else None
    findings = []          # dicts, one per failure
    summary = Counter()    # (source, check, type) -> rows examined
    for src in sources:
        base = load_source(src)
        prot = load_bundle(src, bundle)
        if prot is None:
            continue
        ids = {e for e, typ in prot if typ in ORIGINAL}
        for ext in sorted(ids):
            cpds = aliases.get((src, ext), [])
            # -- row: each protonation row against its own source row
            verdict = {}
            for typ in ORIGINAL:
                k = (ext, typ)
                if k in prot and k in base:
                    summary[(src, "row", typ)] += 1
                    kind, dH, dq = compare(base[k][0], base[k][1], prot[k][0], prot[k][1])
                    verdict[typ] = kind
                    if kind != "ok":
                        findings.append(dict(source=src, external_id=ext, check="row", type=typ,
                                             from_formula=base[k][0], from_charge=base[k][1],
                                             to_formula=prot[k][0], to_charge=prot[k][1],
                                             dH=dH, dq=dq, kind=kind, cpds=cpds))
            for typ in ORIGINAL:            # is the other row a usable fallback?
                other = "SMILE" if typ == "InChI" else "InChI"
                for f in findings:
                    if f["source"] == src and f["external_id"] == ext and f["check"] == "row" and f["type"] == typ:
                        f["other_row"] = verdict.get(other, "absent")
            # -- bundle: InChI row against SMILE row at the same pH
            ki, ks = (ext, "InChI"), (ext, "SMILE")
            if ki in prot and ks in prot:
                summary[(src, "bundle", "InChI~SMILE")] += 1
                kind, dH, dq = compare(prot[ks][0], prot[ks][1], prot[ki][0], prot[ki][1])
                if kind != "ok":
                    findings.append(dict(source=src, external_id=ext, check="bundle", type="InChI~SMILE",
                                         from_formula=prot[ks][0], from_charge=prot[ks][1],
                                         to_formula=prot[ki][0], to_charge=prot[ki][1],
                                         dH=dH, dq=dq, kind=kind, cpds=cpds, other_row=""))
            # -- layers: what the InChI string declares vs what the row stores
            if ki in base:
                summary[(src, "layers", "InChI")] += 1
                lf, lq = inchi_layers(base[ki][2])
                stored_q = as_int(base[ki][1])
                if lf is not None and stored_q is not None and (lq != stored_q or lf != base[ki][0]):
                    findings.append(dict(source=src, external_id=ext, check="layers", type="InChI",
                                         from_formula=lf, from_charge=f"declared {lq:+d}",
                                         to_formula=base[ki][0], to_charge=f"stored {stored_q:+d}",
                                         dH=parse_formula(base[ki][0]).get("H", 0) - parse_formula(lf).get("H", 0),
                                         dq=stored_q - lq,
                                         kind="stored formula/charge disagrees with the InChI's own layers",
                                         cpds=cpds, other_row=""))
            # -- key: the InChIKey row must be the hash of the InChI row
            kk = (ext, "InChIKey")
            if ki in prot and kk in prot:
                summary[(src, "key", "InChIKey")] += 1
                want = inchikey_of(prot[ki][2])
                if want and want != prot[kk][2]:
                    findings.append(dict(source=src, external_id=ext, check="key", type="InChIKey",
                                         from_formula=prot[ki][0], from_charge=prot[ki][1],
                                         to_formula=prot[kk][2], to_charge=want,
                                         dH="", dq="", kind="InChIKey row is not the key of the InChI row",
                                         cpds=cpds, other_row=""))
            # -- source: inchi.tsv against smiles.tsv
            if ki in base and ks in base:
                summary[(src, "source", "InChI~SMILE")] += 1
                kind, dH, dq = compare(base[ks][0], base[ks][1], base[ki][0], base[ki][1])
                if kind != "ok":
                    findings.append(dict(source=src, external_id=ext, check="source", type="InChI~SMILE",
                                         from_formula=base[ks][0], from_charge=base[ks][1],
                                         to_formula=base[ki][0], to_charge=base[ki][1],
                                         dH=dH, dq=dq, kind=kind, cpds=cpds, other_row=""))
    for f in findings:
        f["n_live_reactions"] = sum(rxn.get(c, 0) for c in f["cpds"]) if rxn else ""
    return summary, findings


def report(bundle, summary, findings, allowed):
    print(f"bundle {bundle}")
    print(f"{'source':<9}{'check':<8}{'type':<13}{'rows':>8}{'strict':>7}{'allowed':>9}")
    fails = Counter((f["source"], f["check"], f["type"]) for f in findings if f["kind"] not in INFO_KINDS)
    allow = Counter((f["source"], f["check"], f["type"]) for f in findings
                    if (f["source"], f["external_id"], f["check"]) in allowed)
    for k in sorted(summary):
        print(f"{k[0]:<9}{k[1]:<8}{k[2]:<13}{summary[k]:>8,}{fails[k]:>7}{allow[k]:>9}")
    rowf = [f for f in findings if f["check"] == "row" and f["kind"] not in INFO_KINDS]
    info = [f for f in findings if f["check"] == "row" and f["kind"] in INFO_KINDS]
    if rowf or info:
        print()
        print("row failures by kind (STRICT -- a protonation cannot do this):")
        for kind, n in Counter(f["kind"] for f in rowf).most_common():
            print(f"   {n:>4}  {kind}")
        print("reported, not failed (INFO):")
        for kind, n in Counter(f["kind"] for f in info).most_common():
            print(f"   {n:>4}  {kind}")
        fixable = sum(1 for f in rowf if f.get("other_row") == "ok")
        print(f"\n   {fixable} of {len(rowf)} strict failures have a clean row in the other "
              f"representation; {len(rowf) - fixable} do not.")
        cpds = {c for f in rowf for c in f["cpds"]}
        print(f"   {len(cpds)} ModelSEED compounds are reachable from a strict failure.")


def write_tsv(path, findings):
    cols = ["source", "external_id", "check", "type", "kind", "from_formula", "from_charge",
            "to_formula", "to_charge", "dH", "dq", "other_row", "cpds", "n_live_reactions"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for f in sorted(findings, key=lambda x: (x["check"], x["source"], x["external_id"], x["type"])):
            w.writerow([f.get(c, "") if c != "cpds" else ";".join(f["cpds"]) for c in cols])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", default="marvin_26.1_ph7",
                    help="protonation bundle name without .tsv (default: marvin_26.1_ph7)")
    ap.add_argument("--sources", nargs="*", default=list(SOURCES))
    ap.add_argument("--tsv", help="write one row per failure here")
    ap.add_argument("--reactions", action="store_true",
                    help="count live reactions per failing compound (slower)")
    ap.add_argument("--fail-on-violation", action="store_true",
                    help="exit 1 if any ROW check fails outside the allowlist")
    ap.add_argument("--root", help="validate a Biochemistry/Structures tree at this path instead "
                                   "(e.g. a repaired copy before it is written back)")
    a = ap.parse_args()
    if a.root:
        global STRUCT
        STRUCT = a.root
    allowed = load_allowlist()
    summary, findings = run(a.bundle, a.sources, a.reactions)
    report(a.bundle, summary, findings, allowed)
    if a.tsv:
        write_tsv(a.tsv, findings)
        print(f"\nwrote {len(findings)} findings to {a.tsv}")
    if a.fail_on_violation:
        bad = [f for f in findings if f["check"] == "row" and f["kind"] not in INFO_KINDS
               and (f["source"], f["external_id"], "row") not in allowed]
        if bad:
            print(f"\nFAIL: {len(bad)} protonation rows break dH == dcharge and are not allowlisted",
                  file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
