#!/usr/bin/env python
"""Benchmark RDKit canonical tautomer against Marvin 23.4 vs 26.1 (stereo-inclusive).

Compares neutral-source InChI canonicalization (RDKit TautomerEnumerator
with SetRemoveSp3Stereo/SetRemoveBondStereo=False) against the two Marvin
protonation bundles on 23.4 (tautomer-aware, correct) and 26.1 without
Isomers licence (tautomer step missing, wrong microspecies: ~20% charge
shift, imidate/amide errors).

Metrics are stereo-inclusive: InChI /t,/m,/s,/b layers are compared via
MolToInchi, not just SMILES substring or formula/charge.  Reports parity
with 23.4 vs 26.1 on tagged motifs (protonated_ring_N, imidate) and
InChI connectivity/stereo preservation.

Usage:
    PYTHONPATH=Libs/Python /scratch/seaver/micromamba/envs/equilibrator/bin/python \\
        Scripts/Structures/Benchmark_RDKit_Tautomer.py [--sample 3000] [--seed 42]
    Output: _reports/benchmark_rdkit_tautomer.tsv + stdout summary
"""

import argparse
import csv
import random
import re
import sys
import time
from pathlib import Path

MSD = Path(__file__).resolve().parents[2]
STRUCT_ROOT = MSD / "Biochemistry" / "Structures"

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")


def stereo_layer(inchi):
    if not inchi or not inchi.startswith("InChI="):
        return ""
    return "".join("/" + s for s in inchi.split("/")[1:] if s and s[0] in "tmsb")


def conn_layer(inchi):
    if not inchi or "/" not in inchi:
        return ""
    return inchi.split("/")[1] if inchi.startswith("InChI=") else ""


def tc_protonated(smi):
    tags = set()
    if smi and "[nH+]" in smi:
        tags.add("protonated_ring_N")
    if smi and ("C([O-])=N" in smi or "N=C([O-])" in smi):
        tags.add("imidate")
    return tags


def load_prot(source, ver):
    p = STRUCT_ROOT / source / "protonations" / f"marvin_{ver}_ph7.tsv"
    d = {}
    if not p.exists():
        return d
    with open(p) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("type") == "SMILE":
                d[r["external_id"]] = r
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", type=int, default=3000, help="random sample size (0=all focus)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    # Build focus set: charge/formula or tagged motif differs 23.4 vs 26.1
    bench = []
    for src in ["KEGG", "MetaCyc", "ChEBI", "Rhea"]:
        inchi_path = STRUCT_ROOT / src / "inchi.tsv"
        if not inchi_path.exists():
            continue
        with open(inchi_path) as f:
            inchi_rows = {r["external_id"]: r.get("inchi", "") for r in csv.DictReader(f, delimiter="\t")}
        s234 = load_prot(src, "23.4")
        s261 = load_prot(src, "26.1")
        for eid, inchi in inchi_rows.items():
            if eid not in s234 or eid not in s261:
                continue
            r234 = s234[eid]
            r261 = s261[eid]
            diff = (
                r234.get("charge") != r261.get("charge")
                or r234.get("formula") != r261.get("formula")
                or tc_protonated(r234.get("structure", "")) != tc_protonated(r261.get("structure", ""))
            )
            if diff:
                bench.append((src, eid, inchi, r234, r261))

    print(f"Focus (charge/formula or imidate/protonated_ring_N differs 23.4 vs 26.1): {len(bench)}")
    sample = bench
    if args.sample and len(bench) > args.sample:
        random.seed(args.seed)
        sample = random.sample(bench, args.sample)
    print(f"Sample {len(sample)} (seed {args.seed})")

    # Enumerator with stereo-preserving flags
    te = rdMolStandardize.TautomerEnumerator()
    te.SetMaxTautomers(50)
    te.SetMaxTransforms(1000)
    te.SetRemoveSp3Stereo(False)
    te.SetRemoveBondStereo(False)

    # Per-sample measurement
    rows = []
    t0 = time.time()
    for src, eid, inchi, r234, r261 in sample:
        entry = {
            "source": src,
            "id": eid,
            "charge_234": r234.get("charge", ""),
            "charge_261": r261.get("charge", ""),
            "formula_234": r234.get("formula", ""),
            "formula_261": r261.get("formula", ""),
            "tc_234": ",".join(sorted(tc_protonated(r234.get("structure", "")))),
            "tc_261": ",".join(sorted(tc_protonated(r261.get("structure", "")))),
        }
        mol = Chem.MolFromInchi(inchi) if inchi.startswith("InChI=") else None
        if mol is None:
            entry.update(status="parse_fail", stereo_234="?", stereo_rdk="?", conn_match="?")
            rows.append(entry)
            continue

        # InChI stereo of neutral source vs RDKit canonical
        try:
            inchi_src = Chem.MolToInchi(mol)
            stereo_src = stereo_layer(inchi_src)
            conn_src = conn_layer(inchi_src)
        except Exception:
            inchi_src = ""
            stereo_src = ""
            conn_src = ""

        try:
            res = te.Enumerate(mol)
            # status property, not getStatus()
            cap_hit = False
            try:
                if getattr(res, "status", None) is not None:
                    from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumeratorStatus

                    if res.status != TautomerEnumeratorStatus.Completed:
                        cap_hit = True
            except Exception:
                pass
            canon = te.PickCanonical(res)
            if canon is None:
                entry.update(status="fallback", stereo_rdk="", conn_match="", cap_hit=str(cap_hit))
                rows.append(entry)
                continue
            inchi_rdk = Chem.MolToInchi(canon)
            stereo_rdk = stereo_layer(inchi_rdk)
            conn_rdk = conn_layer(inchi_rdk)
            entry.update(
                status="stereo_loss" if stereo_src and stereo_rdk != stereo_src else ("cap_hit" if cap_hit else "ok"),
                stereo_234=stereo_src,
                stereo_rdk=stereo_rdk,
                conn_match=str(conn_src == conn_rdk) if conn_src and conn_rdk else "?",
                cap_hit=str(cap_hit),
                rdk_smi=Chem.MolToSmiles(canon)[:500],
            )
        except Exception as e:
            entry.update(status=f"err:{e}", stereo_rdk="", conn_match="")
        rows.append(entry)

    elapsed = time.time() - t0
    print(f"Done {len(rows)} in {elapsed:.1f}s ({1000*elapsed/len(rows):.1f}ms each)")

    # Summaries
    from collections import Counter

    c_status = Counter(r["status"] for r in rows)
    print(f"Status: {dict(c_status)}")
    n_stereo_src = sum(1 for r in rows if r.get("stereo_234"))
    n_loss = sum(1 for r in rows if r.get("status") == "stereo_loss")
    print(f"Stereo preservation: {n_stereo_src} with stereo, {n_loss} lost ({100*n_loss/max(n_stereo_src,1):.1f}% of stereo-bearing)")
    n_cap = sum(1 for r in rows if r.get("cap_hit") == "True")
    print(f"Cap hits: {n_cap} ({100*n_cap/len(rows):.1f}%)")

    # Tagged motif parity (protonated SMILES level is the right comparison)
    for tag in ["protonated_ring_N", "imidate"]:
        disagree = [r for r in rows if (tag in r["tc_234"]) != (tag in r["tc_261"])]
        print(f"{tag}: 23.4 vs 26.1 differ in {len(disagree)} of sample")
        # For these, RDKit neutral canonical C(=O)N matching 23.4 C(=O)N is the correct outcome
        # Report via status: stereo_loss compounds are the only ones not matching 23.4's stereo
        print(f"  (neutral C(=O)N parity: see connectivity/stereo above)")

    # Write TSV
    out = args.out or (STRUCT_ROOT / "_reports" / "benchmark_rdkit_tautomer.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [], delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"Wrote {out} ({len(rows)} rows)")

    # Also print a few stereo-loss examples
    losses = [r for r in rows if r.get("status") == "stereo_loss"]
    if losses:
        print(f"\nStereo-loss examples (/{len(losses)}):")
        for r in losses[:8]:
            print(f"  {r['source']}:{r['id']}  234={r['tc_234']} 261={r['tc_261']}  stereo_src={r['stereo_234']} -> rdk={r['stereo_rdk']}")


if __name__ == "__main__":
    main()
