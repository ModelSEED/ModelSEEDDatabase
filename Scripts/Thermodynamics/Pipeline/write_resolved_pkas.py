#!/usr/bin/env python
"""Materialise the resolved pKa table the cache build consumes.

``build_modelseed_cache.py --pkas`` takes one file, but the pKa layer is a
cascade -- Alberty, then IUPAC, then the upstream cache's own ladder, then
MolGpKa -- so the cascade is resolved here
and written out as a single table. Provenance is not lost: the chosen source is
recorded per compound in the companion ``.provenance.tsv``, which is what the
post-build audit reads to prove no ChemAxon-derived value survived.

Kept separate from the builder deliberately. The resolution order is a data
decision that wants inspecting and diffing between runs; the builder is
mechanical.
"""
import os
import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

# Relocated from the eQuilibrator working tree 2026-09-05. ROOT is that
# tree -- it holds the caches and fitted parameters, which are gigabytes and
# are not in this repository -- so it is named by environment variable rather
# than derived from this file's location.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))
sys.path.insert(0, str(Path(__file__).resolve().parent))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "resolved_pkas.tsv")
    ap.add_argument("--cache", type=Path, default=None,
                    help="cache supplying the 'cache' tier's ladders and the "
                         "magnesium guard's reference proton counts. Defaults "
                         "to the upstream Zenodo cache, which is the database "
                         "seed_mapping's our_cache_id indexes; pointing this at "
                         "a cache we rebuilt would be circular.")
    ap.add_argument("--no-guard", action="store_true",
                    help="accept a tier's ladder even if the magnesium guard "
                         "would refuse it, instead of falling through")
    a = ap.parse_args()

    import modelseed_pkas as msp

    structures = msp.load_structures(msp.DEV_DIR)
    ref_cache = a.cache or msp.UPSTREAM_CACHE
    cache_ladders = msp.load_cache_ladders(cache=ref_cache)
    guard = None if a.no_guard else msp.magnesium_guard(ref_cache)
    resolved = msp.resolve_pkas(
        structures,
        msp.load_marvin_pkas(msp.DEV_DIR),
        msp.load_molgpka_pkas(msp.DEV_DIR),
        alberty=msp.load_alberty_pkas(),
        iupac=msp.load_iupac_pkas(),
        cache=cache_ladders,
        admissible=guard,
    )
    print(f"structures {len(structures):,}   resolved {len(resolved):,}   "
          f"cache-tier ladders available {len(cache_ladders):,}"
          f"{'' if guard else '   [guard disabled]'}")
    counts = Counter(src.split(":")[0] for _, src in resolved.values())
    for k, v in counts.most_common():
        print(f"  {k:<10}{v:>8,}  {v / len(resolved):>6.1%}")

    with a.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["external_id", "kind", "value", "tool", "tool_version"])
        for cpd, (values, src) in sorted(resolved.items()):
            w.writerow([cpd, "pKa",
                        ";".join(f"{v:.4f}" for v in values),
                        src.split(":")[0], src])
    prov = a.out.with_suffix(".provenance.tsv")
    with prov.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["seed_id", "source", "detail", "n_pkas"])
        for cpd, (values, src) in sorted(resolved.items()):
            w.writerow([cpd, src.split(":")[0], src, len(values)])
    print(f"\nwrote {a.out}\nwrote {prov}")


if __name__ == "__main__":
    main()
