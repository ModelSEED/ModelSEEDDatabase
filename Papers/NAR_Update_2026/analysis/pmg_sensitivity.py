#!/usr/bin/env python3
"""How far does dG'° move when you change the magnesium condition?

The pipeline predicts at pMg 3.0 (1 mM free Mg2+), inherited rather than
chosen: the pre-2023 ModelSEED scripts set p_h explicitly and never p_mg, so
they silently picked up eQuilibrator's *physiological* default while using the
standard pH. Meanwhile the training loader substitutes pMg 14 -- essentially no
magnesium -- for the 66.4% of TECRDB rows that never recorded it, although the
rows that did record it have a median of 2.4 (about 4 mM).

So the model is fitted mostly as though the assays were magnesium-free, and
then queried at 1 mM. This measures whether that matters, before anyone argues
about which value is right.

Three conditions, same reactions, same parameters:

  pMg 14.0   what two thirds of the training rows are assumed to be
  pMg  3.0   what we currently predict at (eQuilibrator "physiological")
  pMg  2.4   the median of the rows that actually measured it

Reported separately for reactions that touch a Mg2+-binding compound and those
that do not -- the latter are the control and should not move at all.
"""
import os
import argparse, csv, json, glob, sqlite3, sys
from pathlib import Path

# ROOT is the eQuilibrator WORKING TREE, not this repository: these analyses read
# caches and fitted parameters that live there and are far too large to commit.
# It was __file__/../.. while this script lived in that tree; now it must be named
# explicitly. Overridable so the analysis is not pinned to one host.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))
# This repository, derived from the file location rather than hardcoded.
MSD = Path(__file__).resolve().parents[3]
KCAL = 4.184


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=str(ROOT / "data/modelseed_cache_msd/compounds.sqlite"))
    ap.add_argument("--params", default=str(ROOT / "data/cc_params_msd.npz"))
    ap.add_argument("--limit", type=int, default=400)
    ap.add_argument("--out", type=Path, default=ROOT / "data/pmg_sensitivity.tsv")
    a = ap.parse_args()

    import numpy as np
    from component_contribution.parameters import CCModelParameters
    from component_contribution.predict import GibbsEnergyPredictor
    from equilibrator_api import ComponentContribution, Q_
    from equilibrator_cache.api import create_compound_cache_from_sqlite_file

    # which compounds bind Mg2+ at all
    con = sqlite3.connect(f"file:{a.cache}?mode=ro", uri=True)
    seed = con.execute("select id from registries where namespace='seed'").fetchone()[0]
    mg_cids = {c for (c,) in con.execute(
        "select distinct compound_id from magnesium_dissociation_constant")}
    mg_seed = {acc for acc, cid in con.execute(
        "select accession, compound_id from compound_identifiers where registry_id=?",
        (seed,)) if cid in mg_cids}
    con.close()
    print(f"Mg2+-binding compounds with a seed: accession: {len(mg_seed):,}", file=sys.stderr)

    # reactions that already have a usable value, split by Mg involvement
    src = MSD / "Biochemistry/Thermodynamics/eQuilibrator/ModelSEED_Reaction_Energies.tsv"
    with src.open() as fh:
        fh.readline()
        ok = [r for r in csv.DictReader(fh, delimiter="\t") if r["status"] == "ok"]
    stoich = {}
    for f in glob.glob(str(MSD / "Biochemistry/reaction_*.json")):
        for e in json.load(open(f)):
            stoich[e["id"]] = {c["compound"] for c in (e.get("stoichiometry") or [])}
    with_mg = [r for r in ok if stoich.get(r["reaction_id"], set()) & mg_seed]
    without = [r for r in ok if not (stoich.get(r["reaction_id"], set()) & mg_seed)]
    print(f"usable reactions: {len(with_mg):,} touch Mg-binders, {len(without):,} do not",
          file=sys.stderr)
    sample = with_mg[: a.limit] + without[: a.limit]

    cc = ComponentContribution(
        ccache=create_compound_cache_from_sqlite_file(a.cache),
        predictor=GibbsEnergyPredictor(parameters=CCModelParameters.from_npz(Path(a.params))))
    cc.p_h = Q_(7.0); cc.ionic_strength = Q_("0.25 M"); cc.temperature = Q_("298.15 K")

    CONDS = [("pMg14.0", 14.0), ("pMg3.0", 3.0), ("pMg2.4", 2.4)]
    rows = []
    for i, r in enumerate(sample, 1):
        if i % 200 == 0:
            print(f"  {i}/{len(sample)}", file=sys.stderr, flush=True)
        try:
            rxn = cc.parse_reaction_formula(r["formula"])
        except Exception:
            continue
        vals = {}
        for label, pmg in CONDS:
            cc.p_mg = Q_(pmg)
            try:
                vals[label] = cc.standard_dg_prime(rxn).value.magnitude / KCAL
            except Exception:
                vals[label] = None
        if any(v is None for v in vals.values()):
            continue
        rows.append({"reaction_id": r["reaction_id"],
                     "touches_mg": int(bool(stoich.get(r["reaction_id"], set()) & mg_seed)),
                     **{k: f"{v:.4f}" for k, v in vals.items()}})

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    def report(sel, label):
        if not sel:
            return
        d_train = np.array([abs(float(r["pMg3.0"]) - float(r["pMg14.0"])) for r in sel])
        d_meas = np.array([abs(float(r["pMg2.4"]) - float(r["pMg3.0"])) for r in sel])
        print(f"\n{label}  (n={len(sel)})")
        print(f"  |pMg 3.0 - pMg 14|  predict-vs-training-assumption : "
              f"median {np.median(d_train):7.3f}  p90 {np.percentile(d_train,90):8.3f}  "
              f"max {d_train.max():8.3f} kcal/mol")
        print(f"  |pMg 2.4 - pMg 3.0| current-vs-measured-median     : "
              f"median {np.median(d_meas):7.3f}  p90 {np.percentile(d_meas,90):8.3f}  "
              f"max {d_meas.max():8.3f} kcal/mol")
        print(f"  moved >0.5 kcal/mol by the training gap: "
              f"{int((d_train>0.5).sum())} ({100*(d_train>0.5).mean():.0f}%)")

    report([r for r in rows if str(r["touches_mg"]) == "1"], "REACTIONS TOUCHING AN Mg-BINDER")
    report([r for r in rows if str(r["touches_mg"]) == "0"], "CONTROL: no Mg-binder")
    print(f"\nwrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
