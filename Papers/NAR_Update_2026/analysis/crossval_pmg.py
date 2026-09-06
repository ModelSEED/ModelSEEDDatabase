"""Is pMg 14 the right default for the 3,015 rows that never recorded it?

Component contribution's loader substitutes pMg 14 -- essentially zero free
Mg2+ -- for the 66.4% of TECRDB rows with no recorded magnesium. The rows that
DID record it have a median of 2.4, about 4 mM. Enzyme assays for
ATP-dependent chemistry contain magnesium because the chemistry requires it, so
"not recorded" is far more likely to mean "present and unremarkable" than
"absent". Predictions are then made at pMg 3.0 (1 mM), so the model is fitted
mostly as though the assays were magnesium-free and queried as though they were
not.

This has never been tested. The existing `--mg-from-metadata` experiment
(TRAINING_DATA.md, "7e-11 worse, do not use") could not have tested it: it
derives pMg from openTECR's `Cofactor` column, and the base TECRDB has no such
column -- so it only ever moved the ~178 openTECR-sourced rows, inside the
augmented table that was independently harmful. It says nothing about the
default applied to the original 3,015.

Design, mirroring tools/crossval_corrections.py so the numbers are comparable:
folds are drawn once and every arm predicts the SAME held-out measurements;
the only difference is what pMg was assumed for the unrecorded rows during
training. An arm is scored only on reactions well-posed under every arm.

    arm "pmg14"   status quo -- unrecorded rows assumed magnesium-free
    arm "pmg2.4"  unrecorded rows assumed at the median of those that measured
    arm "pmg3.0"  unrecorded rows assumed at the prediction condition, so
                  training and prediction are mutually consistent

A win for 2.4 or 3.0 means the default is not merely inconsistent but wrong.
A null means the transform is insensitive at this scale and the honest fix is
to document the inconsistency rather than re-fit.
"""
import os
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ROOT is the eQuilibrator WORKING TREE, not this repository: these analyses read
# caches and fitted parameters that live there and are far too large to commit.
# It was __file__/../.. while this script lived in that tree; now it must be named
# explicitly. Overridable so the analysis is not pinned to one host.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))
sys.path.insert(0, str(ROOT))
FOLDS = 10
SEED = 20260901


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tecrdb", type=Path,
                    default=ROOT / "data/opentecr/TECRDB_dedup.csv")
    ap.add_argument("--cache", type=Path,
                    default=ROOT / "data/modelseed_cache_msd/compounds.sqlite")
    ap.add_argument("--arms", default="14,2.4,3.0")
    a = ap.parse_args()

    from component_contribution.trainer import Trainer
    from component_contribution.training_data import FullTrainingDataFactory
    from equilibrator_cache import Q_
    from equilibrator_cache.api import create_compound_cache_from_sqlite_file

    # Which rows never recorded pMg? Do NOT try to match rows back to the raw
    # Zenodo table by url -- a url is a NIST page reference shared by many
    # rows, so the ~3,000 unrecorded rows span only 1,207 urls and matching on
    # it also rewrites rows that DID record a value.
    #
    # The filled rows identify themselves: no genuine measurement in the whole
    # corpus exceeds pMg 7.4, so p_mg == 14 exactly means "the loader put this
    # here". Verified against the raw table: 1,529 genuine reports, none at 14.
    src0 = pd.read_csv(a.tecrdb)
    filled = np.isclose(pd.to_numeric(src0["p_mg"], errors="coerce"), 14.0)
    genuine = int((~filled).sum())
    print(f"rows in {a.tecrdb.name}: {len(src0):,}")
    print(f"  pMg genuinely recorded : {genuine:,} "
          f"(median {pd.to_numeric(src0.loc[~filled,'p_mg']).median():g})")
    print(f"  pMg filled with 14     : {int(filled.sum()):,}  <- the arms rewrite these")

    ccache = create_compound_cache_from_sqlite_file(a.cache)
    arms = {}
    for tok in a.arms.split(","):
        val = float(tok)
        src = src0.copy()
        mask = filled
        src.loc[mask, "p_mg"] = val
        p = ROOT / f"data/_pmg_arm_{tok}.csv"
        src.to_csv(p, index=False)
        print(f"  arm pMg {tok:<5} rewrote {int(mask.sum()):,} rows -> {p.name}")
        arms[f"pmg{tok}"] = p

    # Build each arm's training matrices. The stoichiometric matrix is the same
    # across arms (same reactions); only b, the transformed dG, moves.
    mats = {}
    for name, path in arms.items():
        td = FullTrainingDataFactory(ccache=ccache).make(tecrdb_path=path)
        S = td.stoichiometric_matrix
        G = Trainer.group_incidence_matrix(td)
        S.index = S.index.map(lambda c: c.id)
        G.index = G.index.map(lambda c: c.id)
        def mag(series, unit):
            return np.array([float(x.m_as(unit)) if hasattr(x, "m_as") else float(x)
                             for x in series])
        mats[name] = dict(S=S, G=G, b=mag(td.standard_dg, "kJ/mol"),
                          w=mag(td.weight, "dimensionless")
                            if hasattr(td.weight.iloc[0], "m_as")
                            else np.asarray(td.weight, dtype=float),
                          ng=td.group_summary.shape[0])
        print(f"  {name}: S {S.shape}, b range "
              f"{mats[name]['b'].min():.1f} .. {mats[name]['b'].max():.1f} kJ/mol")

    ref = mats[next(iter(mats))]
    ncols = ref["S"].shape[1]
    if any(m["S"].shape[1] != ncols for m in mats.values()):
        sys.exit("FATAL: arms have different reaction counts -- the test set "
                 "would not be identical and the comparison is meaningless")

    # How much did the assumption actually move the training targets?
    base = mats["pmg14.0"]["b"] if "pmg14.0" in mats else mats["pmg14"]["b"]
    for name, m in mats.items():
        d = np.abs(m["b"] - base)
        print(f"  |b - b(pMg14)| for {name:<8} median {np.median(d):7.4f}  "
              f"max {d.max():8.4f} kJ/mol   moved>0.1: {int((d>0.1).sum()):,}")

    rng = np.random.default_rng(SEED)
    folds = np.array_split(rng.permutation(np.arange(ncols)), FOLDS)

    def fit_predict(m, train_cols, test_cols):
        S, G, b, w = m["S"], m["G"], m["b"], m["w"]
        Ss = S.iloc[:, train_cols]
        keep = (Ss != 0).any(axis=1)
        Ss, Gs = Ss.loc[keep], G.loc[keep]
        p = Trainer.train_from_matrices(
            Ss, Gs,
            pd.Series([Q_(float(x), "kJ/mol") for x in b[train_cols]], index=Ss.columns),
            pd.Series(w[train_cols], index=Ss.columns), m["ng"])
        dgf = pd.Series(np.asarray(p.dG0_cc).flatten(), index=Ss.index)
        return S.iloc[:, test_cols].loc[keep].T.values @ dgf.values - b[test_cols]

    def wellposed(m, train_cols, test_cols):
        present = (m["S"].iloc[:, train_cols] != 0).any(axis=1)
        return np.array([c for c in test_cols
                         if bool(present[m["S"].iloc[:, c] != 0].all())], dtype=int)

    res = {k: [] for k in mats}
    for i, fold in enumerate(folds, 1):
        train = np.setdiff1d(np.arange(ncols), fold)
        scored = fold
        for m in mats.values():
            scored = np.intersect1d(scored, wellposed(m, train, fold))
        if not len(scored):
            continue
        line = f"  fold {i}/{FOLDS}: held out {len(scored):4d} "
        for name, m in mats.items():
            r = fit_predict(m, train, scored)
            res[name].append(r)
            line += f"  {name} {np.sqrt((r**2).mean()):6.3f}"
        print(line, flush=True)

    print("\ncross-validated prediction of held-out measurements "
          "(identical test set in every arm):")
    out = {}
    for name in mats:
        r = np.concatenate(res[name])
        out[name] = r
        print(f"  {name:<9} RMSE {np.sqrt((r**2).mean()):7.4f}   "
              f"MAE {np.abs(r).mean():7.4f}   median|e| {np.median(np.abs(r)):7.4f}")
    try:
        from scipy.stats import wilcoxon
        base_name = "pmg14.0" if "pmg14.0" in out else "pmg14"
        for name in out:
            if name == base_name:
                continue
            s, p = wilcoxon(np.abs(out[base_name]), np.abs(out[name]))
            better = np.median(np.abs(out[name])) < np.median(np.abs(out[base_name]))
            verdict = ("no significant difference" if p >= 0.05
                       else ("BETTER than pMg 14" if better else "WORSE than pMg 14"))
            print(f"  {name:<9} vs {base_name}: p = {p:.4g}   {verdict}")
    except ImportError:
        print("  scipy unavailable; no significance test")


if __name__ == "__main__":
    main()
