"""Generate formation energies for every ModelSEED compound.

Posed exactly as ModelSEED's 2020 pipeline posed it
(Scripts/Thermodynamics/Retrieve_eQuilibrator_Compound_Energies.py): the
formation energy is standard_dg_prime of the half reaction " = <compound>",
at pH 7 / ionic strength 0.25 M / pMg 3.0 / 298.15 K.

Output carries both kJ/mol and kcal/mol. ModelSEED stores kcal/mol.
"""

import os
import argparse
import csv
import glob
import sys
import warnings
from collections import Counter
from pathlib import Path

warnings.simplefilter("ignore")
# Relocated from the eQuilibrator working tree 2026-09-05. ROOT is that
# tree -- it holds the caches and fitted parameters, which are gigabytes and
# are not in this repository -- so it is named by environment variable rather
# than derived from this file's location.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
MSD = Path("/scratch/seaver/Claude_Projects/MSD_Structures/ModelSEEDDatabase/Biochemistry")
KCAL = 4.184


def load_compounds():
    out = []
    for f in sorted(glob.glob(str(MSD / "compound_*.tsv"))):
        if ".provenance." in f:
            continue
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r.get("id"):
                out.append((r["id"], (r.get("name") or "").strip(),
                            (r.get("formula") or "").strip()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path,
                    default=ROOT / "data/modelseed_cache_b/compounds.sqlite")
    ap.add_argument("--params", default=str(ROOT / "data/cc_params_dedup.npz"))
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data/modelseed_formation_energies.tsv")
    ap.add_argument("--p-h", type=float, default=7.0)
    ap.add_argument("--ionic-strength", type=float, default=0.25)
    ap.add_argument("--p-mg", type=float, default=3.0)
    ap.add_argument("--temperature", type=float, default=298.15)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    from component_contribution.parameters import CCModelParameters
    from component_contribution.predict import GibbsEnergyPredictor
    # equilibrator_api.Reaction subclasses the cache's Reaction and adds the
    # methods standard_dg_prime needs (separate_stored_dg_prime). Importing
    # the cache class instead fails with an AttributeError at prediction time.
    from equilibrator_api import ComponentContribution, Q_, Reaction
    from equilibrator_cache.api import create_compound_cache_from_sqlite_file

    kw = {}
    if str(args.params).lower() != "none":
        kw["predictor"] = GibbsEnergyPredictor(
            parameters=CCModelParameters.from_npz(Path(args.params)))
    cc = ComponentContribution(
        ccache=create_compound_cache_from_sqlite_file(args.cache), **kw)
    cc.p_h = Q_(args.p_h)
    cc.ionic_strength = Q_(f"{args.ionic_strength} M")
    cc.p_mg = Q_(args.p_mg)
    cc.temperature = Q_(f"{args.temperature} K")

    RMSE_INF = float(cc.predictor.preprocess.RMSE_inf)
    cpds = load_compounds()
    if args.limit:
        cpds = cpds[:args.limit]
    print(f"ModelSEED compounds: {len(cpds)}", flush=True)

    stats = Counter()
    with args.out.open("w", newline="") as fh:
        fh.write(f"# cache={args.cache} params={args.params} "
                 f"p_h={args.p_h} ionic_strength={args.ionic_strength}M "
                 f"p_mg={args.p_mg} T={args.temperature}K\n")
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["compound_id", "name", "formula", "status",
                    "dgf_prime_kJ_per_mol", "uncertainty_kJ_per_mol",
                    "dgf_prime_kcal_per_mol", "uncertainty_kcal_per_mol"])
        for i, (cid, name, formula) in enumerate(cpds, 1):
            if i % 5000 == 0:
                print(f"  {i}/{len(cpds)}", flush=True)
            try:
                rxn = Reaction.parse_formula(cc.get_compound, f" = seed:{cid}")
            except Exception:
                stats["compound not in cache"] += 1
                w.writerow([cid, name, formula, "compound not in cache",
                            "", "", "", ""])
                continue
            # Same two-branch split as the reaction generator -- see the long
            # comment there. `residual` non-empty means the compound cannot be
            # expressed in the model's basis at all: mu is discarded and the
            # return is literally zero plus the pH transform. Anything else is
            # a real estimate, however wide its error bar, and is published.
            # Split off compounds carrying a measured dG before asking, exactly
            # as ComponentContribution.standard_dg_prime does internally --
            # otherwise a measured-but-undecomposable compound is reported as
            # undecomposable and we decline a value the library would give.
            try:
                _residual_rxn, _stored = rxn.separate_stored_dg_prime(
                    p_h=cc.p_h, p_mg=cc.p_mg,
                    ionic_strength=cc.ionic_strength,
                    temperature=cc.temperature)
                _mu, _s_fin, _s_inf, residual = \
                    cc.predictor.get_reaction_prediction(_residual_rxn)
            except Exception:
                stats["prediction failed"] += 1
                w.writerow([cid, name, formula, "prediction failed", "", "", "", ""])
                continue

            if residual:
                stats["undecomposable"] += 1
                w.writerow([cid, name, formula, "undecomposable", "", "", "", ""])
                continue

            try:
                dg = cc.standard_dg_prime(rxn)
            except Exception:
                stats["prediction failed"] += 1
                w.writerow([cid, name, formula, "prediction failed", "", "", "", ""])
                continue

            mu, err = dg.value.magnitude, dg.error.magnitude
            stats["ok"] += 1
            # True sigma, never the string "inf" -- see the reaction generator.
            if err >= 0.5 * RMSE_INF:
                stats["  ...of which effectively unconstrained"] += 1
            elif err > 100.0:
                stats["  ...of which very low confidence"] += 1
            w.writerow([cid, name, formula, "ok", f"{mu:.3f}", f"{err:.3f}",
                        f"{mu / KCAL:.3f}", f"{err / KCAL:.3f}"])

    print(f"\nwrote {args.out}")
    tot = sum(stats.values())
    for k, v in stats.most_common():
        print(f"  {k:<26} {v:6d}  ({100*v/tot:4.1f}%)")


if __name__ == "__main__":
    main()
