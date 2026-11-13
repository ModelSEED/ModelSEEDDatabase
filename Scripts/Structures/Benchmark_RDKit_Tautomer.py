#!/usr/bin/env python
"""Benchmark the shipped RDKit tautomer step (iii-RDKit) against Marvin 23.4.

Focus set: compounds where the 26.1 bundle generated WITHOUT step (iii)
disagrees with 23.4 on charge/formula or on the imidate/protonated-ring-N
motifs.  For each, runs the SAME tautomer_canonical() that
Run_Marvin_Protonations.py ships (uncharged first, capped canonical kept,
stereo guard) on the source InChI, and measures parity of the structure
handed to Marvin with 23.4's shipped protonated SMILES:

  - conn_match_src / conn_match_234: InChI connectivity of the tautomer
    vs the source, and vs 23.4's protonated output
  - motif parity: amide vs imidate of the tautomer vs 23.4's protonated
    motif; the imidate target table (fixed / still imidate) is the
    headline metric for the defect step (iii) was restored to fix
  - stereo: /t,/m,/s,/b layers of source (stereo_src), of 23.4's output
    (stereo_234), and of the tautomer (stereo_rdk); losses are counted
    and reported both as a share of all rows and of stereo-bearing rows

This does NOT run Marvin 26.1 (no Isomers/Protonation licence here): it
compares the input handed to MajorMicrospeciesPlugin, which performs no
tautomer shifts, so an amide handed in stays an amide.  Final-charge
parity against 23.4 requires a licensed run of Run_Marvin_Protonations.py
over the sample; nothing here claims otherwise.

Usage:
    PYTHONPATH=Libs/Python /scratch/seaver/micromamba/envs/equilibrator/bin/python \
        Scripts/Structures/Benchmark_RDKit_Tautomer.py [--sample 3000] [--seed 42]
    Output: _reports/benchmark_rdkit_tautomer.tsv + stdout summary
"""

import argparse
import csv
import os
import random
import sys
import time
from collections import Counter
from pathlib import Path

MSD = Path(__file__).resolve().parents[2]
STRUCT_ROOT = MSD / "Biochemistry" / "Structures"
sys.path.insert(0, str(MSD / "Scripts" / "Structures"))

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

from Run_Marvin_Protonations import tautomer_canonical

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
    """Motif tags of a PROTONATED SMILES (23.4/26.1 bundle rows)."""
    tags = set()
    if smi and "[nH+]" in smi:
        tags.add("protonated_ring_N")
    if smi and ("C([O-])=N" in smi or "N=C([O-])" in smi):
        tags.add("imidate")
    return tags


def tc_neutral(smi):
    """Motif tags of a NEUTRAL tautomer SMILES (what RDKit hands Marvin)."""
    tags = set()
    if smi and ("C(=O)N" in smi or "C(N)=O" in smi or "NC(=O)" in smi):
        tags.add("amide")
    if smi and ("N=C(O)" in smi or "C(O)=N" in smi):
        tags.add("imidic_acid")
    if smi and ("N=C([O-])" in smi or "C([O-])=N" in smi):
        tags.add("imidate_anion")
    if smi and "[nH+]" in smi:
        tags.add("protonated_ring_N")
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

    # Same enumerator/uncharger configuration as Protonator.__init__
    te = rdMolStandardize.TautomerEnumerator()
    te.SetMaxTautomers(50)
    te.SetMaxTransforms(1000)
    te.SetRemoveSp3Stereo(False)
    te.SetRemoveBondStereo(False)
    uncharger = rdMolStandardize.Uncharger()

    rows = []
    t0 = time.time()
    for i, (src, eid, inchi, r234, r261) in enumerate(sample):
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
            entry.update(status="parse_fail", cap_hit="False", stereo_src="?",
                         stereo_234="?", stereo_rdk="?", conn_match_src="?",
                         conn_match_234="?", motif_rdk="?", rdk_smi="")
            rows.append(entry)
            continue

        try:
            inchi_src = Chem.MolToInchi(mol)
            stereo_src = stereo_layer(inchi_src)
            conn_src = conn_layer(inchi_src)
        except Exception:
            inchi_src, stereo_src, conn_src = "", "", ""

        # 23.4's shipped protonated output: its stereo and connectivity
        p234 = r234.get("structure", "")
        pm = Chem.MolFromSmiles(p234) if p234 else None
        try:
            stereo_234 = stereo_layer(Chem.MolToInchi(pm)) if pm is not None else ""
            conn_234 = conn_layer(Chem.MolToInchi(pm)) if pm is not None else ""
        except Exception:
            stereo_234, conn_234 = "", ""

        # The shipped step: uncharge -> enumerate -> capped canonical kept
        out_mol, info = tautomer_canonical(mol, te, uncharger, Chem)
        cap_hit = info == "tautomer_cap"
        try:
            inchi_rdk = Chem.MolToInchi(out_mol) if out_mol is not None else ""
            stereo_rdk = stereo_layer(inchi_rdk)
            conn_rdk = conn_layer(inchi_rdk)
            rdk_smi = Chem.MolToSmiles(out_mol)[:500] if out_mol is not None else ""
        except Exception:
            stereo_rdk, conn_rdk, rdk_smi = "", "", ""
        entry.update(
            status=info,
            cap_hit=str(cap_hit),
            stereo_src=stereo_src,
            stereo_234=stereo_234,
            stereo_rdk=stereo_rdk,
            conn_match_src=str(conn_src == conn_rdk) if conn_src and conn_rdk else "?",
            conn_match_234=str(conn_rdk == conn_234) if conn_rdk and conn_234 else "?",
            motif_rdk=",".join(sorted(tc_neutral(rdk_smi))),
            rdk_smi=rdk_smi,
        )
        rows.append(entry)
        if i % 500 == 0:
            print(f"  {i}/{len(sample)} ({time.time()-t0:.0f}s)", flush=True)

    elapsed = time.time() - t0
    print(f"Done {len(rows)} in {elapsed:.1f}s ({1000*elapsed/len(rows):.1f}ms each)")

    # Summaries
    c_status = Counter(r["status"] for r in rows)
    print(f"Status: {dict(c_status)}")
    n_cap = sum(1 for r in rows if r.get("cap_hit") == "True")
    print(f"Cap hits (counted, canonical still used): {n_cap} ({100*n_cap/len(rows):.1f}%)")
    n_stereo_src = sum(1 for r in rows if r.get("stereo_src"))
    n_loss = sum(1 for r in rows if r.get("status") == "stereo_loss")
    print(f"Stereo loss: {n_loss}/{len(rows)} = {100*n_loss/len(rows):.1f}% of all rows, "
          f"{100*n_loss/max(n_stereo_src,1):.1f}% of the {n_stereo_src} stereo-bearing rows")
    for col in ("conn_match_src", "conn_match_234"):
        c = Counter(r[col] for r in rows)
        print(f"{col}: {dict(c)}")

    # Headline parity vs 23.4: the imidate defect this step exists to fix
    targets = [r for r in rows
               if "imidate" in r["tc_261"].split(",") and "imidate" not in r["tc_234"].split(",")]
    fixed = [r for r in targets
             if "amide" in r["motif_rdk"].split(",") and "imidate_anion" not in r["motif_rdk"].split(",")]
    still = len(targets) - len(fixed)
    print(f"Imidate targets (26.1 wrote imidate, 23.4 wrote amide): {len(targets)}")
    print(f"  fixed by (iii-RDKit): {len(fixed)} ({100*len(fixed)/max(len(targets),1):.1f}%)")
    print(f"  still imidate:        {still} ({100*still/max(len(targets),1):.1f}%)")

    # Write TSV
    out = args.out or (STRUCT_ROOT / "_reports" / "benchmark_rdkit_tautomer.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [], delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"Wrote {out} ({len(rows)} rows)")

    losses = [r for r in rows if r.get("status") == "stereo_loss"]
    if losses:
        print(f"\nStereo-loss examples ({len(losses)}; the UNCHARGED INPUT is handed to Marvin,")
        print("so stereo_rdk below equals stereo_src by construction):")
        for r in losses[:8]:
            print(f"  {r['source']}:{r['id']}  stereo_src={r['stereo_src'][:60]}")


if __name__ == "__main__":
    main()
