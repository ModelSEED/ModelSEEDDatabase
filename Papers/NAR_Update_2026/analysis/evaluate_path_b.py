"""Evaluate the original / Path A / Path B caches over real ModelSEED reactions.

Answers the practical question: for the reactions ModelSEED actually contains,
how many get a usable Delta-G, how uncertain is it, and does anchoring on our
own structures change the answer?
"""

import os
import argparse
import csv
import glob
import re
import sys
import warnings
from collections import Counter
from pathlib import Path

warnings.simplefilter("ignore")
# ROOT is the eQuilibrator WORKING TREE, not this repository: these analyses read
# caches and fitted parameters that live there and are far too large to commit.
# It was __file__/../.. while this script lived in that tree; now it must be named
# explicitly. Overridable so the analysis is not pinned to one host.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))
sys.path.insert(0, str(ROOT))

MSD = str(Path(__file__).resolve().parents[3] / "Biochemistry")


def load_reactions(limit=None, priority_only=False):
    """Parse ModelSEED stoichiometry into seed: reaction formulas."""
    priority = set()
    if priority_only:
        p = Path(__file__).resolve().parents[4] / "priority" / "template_reactions.txt"
        if p.is_file():
            priority = {l.strip() for l in p.open() if l.strip()}

    out = []
    for f in sorted(glob.glob(f"{MSD}/reaction_*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            rid = r.get("id")
            if not rid or (priority and rid not in priority):
                continue
            if str(r.get("is_transport", "0")).strip() in ("1", "true", "True"):
                continue
            stoich = r.get("stoichiometry") or ""
            left, right = [], []
            ok = True
            for term in stoich.split(";"):
                parts = term.split(":")
                if len(parts) < 3:
                    ok = False
                    break
                try:
                    coeff = float(parts[0])
                except ValueError:
                    ok = False
                    break
                cpd, comp = parts[1], parts[2]
                if comp != "0":          # cytosol only; skip multi-compartment
                    ok = False
                    break
                side = right if coeff > 0 else left
                n = abs(coeff)
                side.append(f"{n:g} seed:{cpd}" if n != 1 else f"seed:{cpd}")
            if ok and left and right:
                out.append((rid, " + ".join(left) + " = " + " + ".join(right)))
            if limit and len(out) >= limit:
                return out
    return out


def evaluate(name, cache_path, params_path, reactions):
    from component_contribution.parameters import CCModelParameters
    from component_contribution.predict import GibbsEnergyPredictor
    from equilibrator_api import ComponentContribution
    from equilibrator_cache.api import create_compound_cache_from_sqlite_file

    kw = {}
    if params_path:
        kw["predictor"] = GibbsEnergyPredictor(
            parameters=CCModelParameters.from_npz(params_path))
    cc = ComponentContribution(
        ccache=create_compound_cache_from_sqlite_file(cache_path), **kw)

    res, stats = {}, Counter()
    for rid, formula in reactions:
        try:
            rxn = cc.parse_reaction_formula(formula)
        except Exception:
            stats["unparseable"] += 1
            continue
        if None in rxn.sparse:
            stats["compound missing"] += 1
            continue
        try:
            dg = cc.standard_dg_prime(rxn)
        except Exception:
            stats["prediction failed"] += 1
            continue
        mu, err = dg.value.magnitude, dg.error.magnitude
        if err > 1e4:
            stats["outside CC span (infinite uncertainty)"] += 1
        else:
            stats["usable"] += 1
            res[rid] = (mu, err)
    return name, res, stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=3000)
    ap.add_argument("--priority-only", action="store_true")
    args = ap.parse_args()

    reactions = load_reactions(limit=args.limit, priority_only=args.priority_only)
    print(f"evaluating {len(reactions)} ModelSEED reactions\n")

    configs = [
        ("original", Path.home() / ".cache/equilibrator/compounds.sqlite", None),
        ("path A", ROOT / "data/modelseed_cache/compounds.sqlite", None),
        ("path B", ROOT / "data/modelseed_cache_b/compounds.sqlite",
         ROOT / "data/cc_params_stage2.npz"),
        # Same cache as path B, but trained on the de-duplicated TECRDB plus
        # the 144 structure-first openTECR rows.
        ("dedup", ROOT / "data/modelseed_cache_b/compounds.sqlite",
         ROOT / "data/cc_params_dedup.npz"),
        ("best", ROOT / "data/modelseed_cache_b/compounds.sqlite",
         ROOT / "data/cc_params_best.npz"),
    ]
    results = {}
    for name, cache, params in configs:
        n, res, stats = evaluate(name, cache, params, reactions)
        results[n] = res
        total = sum(stats.values())
        usable = stats["usable"]
        med = sorted(e for _m, e in res.values())
        print(f"{n:<10} usable {usable:5d}/{total:<5d} ({100*usable/max(total,1):4.1f}%)"
              f"   median sigma {med[len(med)//2]:6.2f} kJ/mol" if med else f"{n}: none")
        for k, v in stats.most_common():
            if k != "usable":
                print(f"             {k:<40} {v}")

    # Median sigma across configs is confounded by coverage: a config that
    # rescues additional reactions is scored on a harder set, so its median
    # can rise while no individual estimate got worse. Compare like for like
    # on the reactions every config can compute.
    common = set.intersection(*(set(r) for r in results.values())) if results else set()
    print(f"\nlike-for-like on the {len(common)} reactions usable in ALL configs:")
    for n, res in results.items():
        errs = sorted(res[r][1] for r in common)
        if not errs:
            continue
        med = errs[len(errs) // 2]
        mean = sum(errs) / len(errs)
        print(f"   {n:<10} median sigma {med:6.3f}   mean {mean:6.3f} kJ/mol")

    a, b = results.get("path A", {}), results.get("path B", {})
    both = set(a) & set(b)
    diffs = [(abs(a[r][0] - b[r][0]), r) for r in both]
    changed = [d for d in diffs if d[0] > 0.5]
    print(f"\npath A vs path B over {len(both)} shared reactions:")
    print(f"   |delta| > 0.5 kJ/mol : {len(changed)} ({100*len(changed)/max(len(both),1):.1f}%)")
    if changed:
        changed.sort(reverse=True)
        print(f"   largest deltas: " +
              ", ".join(f"{r}:{d:.1f}" for d, r in changed[:6]))
    tighter = sum(1 for r in both if b[r][1] < a[r][1] - 1e-9)
    looser = sum(1 for r in both if b[r][1] > a[r][1] + 1e-9)
    print(f"   uncertainty tighter under B: {tighter}   looser: {looser}")
    only_b = set(b) - set(a)
    only_a = set(a) - set(b)
    print(f"   usable only under B: {len(only_b)}   only under A: {len(only_a)}")


if __name__ == "__main__":
    main()
