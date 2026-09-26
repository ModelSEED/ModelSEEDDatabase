#!/usr/bin/env python
"""#35 -- which compounds genuinely CANNOT be created, and so must keep a
Zenodo carry-over row?

The plan for a fully open-source pKa layer replaces the 22,326 carried-over
compounds by building them from ModelSEED structures. Carry-over is
load-bearing rather than an optimisation, though: eQuilibrator's group
decomposer rejects small inorganics outright -- phosphate, water, CO2, ammonia
-- and several of those are training anchors. Until now the size of that set
has been an assumption ("tens, not thousands"). This measures it.

A compound must stay carried over if any of:
  * no ModelSEED structure at all
  * the group decomposer refuses the structure
  * it is a training-only foreign accession (metanetx / coco / unresolved kegg),
    which is inherent -- the measurements are KEGG-keyed, not ours
"""
import sys
import os
import argparse, csv, sys, time
from pathlib import Path

# Relocated from the eQuilibrator working tree 2026-09-05. ROOT is that
# tree -- it holds the caches and fitted parameters, which are gigabytes and
# are not in this repository -- so it is named by environment variable rather
# than derived from this file's location.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
MAPPING = ROOT / "data" / "seed_mapping.tsv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "must_carry_over.tsv")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    import modelseed_pkas as msp

    rows = list(csv.DictReader(MAPPING.open(), delimiter="\t"))
    carry = [r for r in rows if (r["our_cache_id"] or "").strip()]
    print(f"carried over today: {len(carry):,}", flush=True)
    if a.limit:
        carry = carry[: a.limit]

    out, t0 = [], time.time()
    counts = {"no structure": 0, "decomposer refuses": 0, "creatable": 0}
    for i, r in enumerate(carry, 1):
        smi = (r["our_smiles"] or "").strip()
        if not smi:
            counts["no structure"] += 1
            out.append((r["seed_id"], "no_structure", "", r["formula"]))
        elif not msp.is_group_decomposable(smi):
            counts["decomposer refuses"] += 1
            out.append((r["seed_id"], "not_decomposable", smi[:80], r["formula"]))
        else:
            counts["creatable"] += 1
        if i % 2000 == 0:
            el = time.time() - t0
            print(f"  {i:,}/{len(carry):,}  {el/i*1000:.0f} ms/cpd  "
                  f"ETA {(len(carry)-i)*el/i/60:.0f} min", flush=True)

    with a.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["seed_id", "reason", "smiles", "formula"])
        w.writerows(out)
    print(f"\nwrote {a.out}  ({len(out):,} rows)")
    for k, v in counts.items():
        print(f"  {k:<22}{v:>8,}")
    print(f"\n=> {counts['creatable']:,} of {len(carry):,} carry-overs "
          f"({counts['creatable']/len(carry):.1%}) could be built from our own "
          f"structures instead.")


if __name__ == "__main__":
    main()
