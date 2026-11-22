#!/usr/bin/env python
"""Heterocycle-N / corrin-porphyrin override: Marvin 26.1 protonates ring
nitrogens at pH 7 that 23.4 (and measurement) leave neutral.

SCOPE. Marvin 26.1 assigns basic pKa ~6.0-8.0 to aromatic ring N in
imidazole/purine/benzimidazole heterocycles and to the central N of corrin
and porphyrin macrocycles, so its pH-7 major microspecies carries [nH+]/[N+]
where 23.4 does not. Measured counterexample: L-histidine imidazole N,
IUPAC pKaH 6.04 (23.4: 6.14, 26.1: 7.20); histamine 6.05 (23.4: 6.14,
26.1: 7.20). The adenine N1 case is handled separately by Adenine_N1.py and
is NOT touched here.

Only structures whose 26.1 bundle charge differs from the 23.4 bundle charge
AND whose 26.1 ladder has a basic site in pKb 6.0-8.0 absent from 23.4's
ladder are processed (the audit list, 134 source ids). Within those, only
ring N atoms carrying +1 are deprotonated: with an H (remove the H, set
charge 0, adenine-style) or without (corrin/porphyrin central N coordinated
to metal: set charge 0 only -- removing an H would be wrong, there is none).

Usage: Heterocycle_N_Override.py [--dry-run] [--sources KEGG MetaCyc ...]
Writes corrected rows to _reports/ and prints counts. Does NOT rewrite the
bundles in place; application is a separate step.
"""

import csv
import os
import sys
from collections import defaultdict

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

STRUCT_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "..", "Biochemistry", "Structures")
SOURCES = ["ChEBI", "KEGG", "MetaCyc", "Rhea"]
AUDIT = "/tmp/near7_261only.tsv"   # set by caller; ext ids with 26.1-only near-7 pKb
PKb_LO, PKb_HI = 6.0, 8.0


def load_bundle_smiles(src, ver):
    d = {}
    with open(os.path.join(STRUCT_ROOT, src, "protonations",
                           f"marvin_{ver}_ph7.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["type"] == "SMILE":
                d[r["external_id"]] = r
    return d


def load_pkb(src, ver):
    d = defaultdict(list)
    path = os.path.join(STRUCT_ROOT, src, "pkas", f"marvin_{ver}.tsv")
    if not os.path.exists(path):
        return d
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["kind"] != "pKb":
                continue
            for tok in (r.get("value") or "").split(";"):
                p = tok.split(":")
                if len(p) == 3:
                    try:
                        d[r["external_id"]].append(float(p[2]))
                    except ValueError:
                        pass
    return d


def audit_targets(sources):
    """ext ids (per source) whose 26.1 charge != 23.4 charge and which have a
    26.1-only basic site in [PKb_LO, PKb_HI]."""
    targets = defaultdict(set)
    for src in sources:
        b4 = load_bundle_smiles(src, "23.4")
        b6 = load_bundle_smiles(src, "26.1")
        p4 = load_pkb(src, "23.4")
        p6 = load_pkb(src, "26.1")
        for eid in set(b4) & set(b6):
            if b4[eid]["charge"] == b6[eid]["charge"]:
                continue
            # only the direction 26.1 adds a proton (more positive charge)
            try:
                c4, c6 = int(b4[eid]["charge"]), int(b6[eid]["charge"])
            except (ValueError, TypeError):
                continue
            if c6 <= c4:
                continue
            kb4 = [v for v in p4.get(eid, []) if PKb_LO <= v <= PKb_HI]
            kb6 = [v for v in p6.get(eid, []) if PKb_LO <= v <= PKb_HI]
            if kb6 and not kb4:
                targets[src].add(eid)
    return targets


def deprotonate_ring_n(smiles, max_sites=1):
    """Deprotonate up to max_sites ring N atoms carrying +1 (the ones 26.1
    added vs 23.4), chosen as the ring N+ with the most N+ neighbours (the
    macrocycle central N in corrins/porphyrins). Returns (new_smiles, n)."""
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    if m is None:
        return smiles, 0
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    rw = Chem.RWMol(m)
    cand = [a.GetIdx() for a in rw.GetAtoms()
            if a.GetSymbol() == "N" and a.GetFormalCharge() == 1 and a.IsInRing()]
    if not cand:
        return smiles, 0
    def nplus_neighbors(i):
        a = rw.GetAtomWithIdx(i)
        return sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "N" and n.GetFormalCharge() == 1)
    cand.sort(key=lambda i: (-nplus_neighbors(i), i))
    targets = cand[:max_sites]
    for i in sorted(targets, reverse=True):
        a = rw.GetAtomWithIdx(i)
        if a.GetTotalNumHs(includeNeighbors=True) >= 1:
            a.SetNumExplicitHs(max(a.GetTotalNumHs() - 1, 0))
        a.SetFormalCharge(0)
    return Chem.MolToSmiles(rw.GetMol()), len(targets)


def main():
    dry = "--dry-run" in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    sources = args or SOURCES
    targets = audit_targets(sources)
    total_cpds = sum(len(v) for v in targets.values())
    print(f"audit targets: {total_cpds} compounds across {dict((k, len(v)) for k, v in targets.items())}")
    grand_sites = 0
    report = []
    for src in sources:
        b6 = load_bundle_smiles(src, "26.1")
        b4 = load_bundle_smiles(src, "23.4")
        for eid in sorted(targets[src]):
            row = b6.get(eid)
            if not row:
                continue
            try:
                diff = int(row["charge"]) - int(b4[eid]["charge"])
            except (ValueError, TypeError):
                continue
            new_smi, n = deprotonate_ring_n(row["structure"], max_sites=max(diff, 1))
            if n:
                grand_sites += n
                report.append((src, eid, row["charge"], new_smi, n))
    print(f"structures to override: {len(report)}, ring N sites deprotonated: {grand_sites}")
    if not dry:
        rep = os.path.join(STRUCT_ROOT, "_reports", "heterocycle_n_override.tsv")
        with open(rep, "w") as f:
            f.write("source\texternal_id\tcharge_261\tnew_smiles\tn_sites\n")
            for r in report:
                f.write("\t".join([r[0], r[1], r[2], r[3], str(r[4])]) + "\n")
        print(f"wrote {rep}")


if __name__ == "__main__":
    main()
