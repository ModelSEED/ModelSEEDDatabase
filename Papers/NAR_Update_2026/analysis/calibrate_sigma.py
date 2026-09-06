"""What does each source's sigma actually mean? One empirical scale per source.

A sigma of 10 kcal/mol does not mean the same thing in the three columns of
the thermodynamics dict:

  Group contribution  propagated uncertainty on the Jankowski 2008 group
                      energies, summed in quadrature over a molecule's groups
  eQuilibrator        fitted covariance of the component-contribution model,
                      ||sigma_fin|| + RMSE_inf*||sigma_inf||, with the second
                      term acting as an explicit "unconstrained" sentinel
  dGPredictor         BayesianRidge posterior predictive standard deviation

They are not interchangeable and cannot be put on a common axis by rescaling:
an earlier attempt on GC showed the miscalibration is a function of sigma
itself (factor 1.17 in the smallest quartile, 0.27 in the largest) and
non-Gaussian, so no affine transform makes the distributions comparable.

What CAN be done is give each source its own empirically measured scale: for
sigma in band X, the observed error is Y. That needs held-out predictions with
their predicted sigma, which is what this produces.

  GC   measured directly against TECRDB. GC never trains on it, so every
       comparison is already out-of-sample -- no folding needed.
  eQ   10-fold over the component-contribution training set. Each fold refits
       from the training columns only, then asks the LIBRARY's own preprocessor
       for mu and sigma on the held-out reactions, so the sigma is the one
       eQuilibrator would really report, not a reimplementation.
  dG   10-fold over the dGPredictor design matrix with BayesianRidge's
       return_std, which is the same quantity the shipped model reports.

Emits data/sigma_calibration_{eq,dg}.tsv with one row per held-out prediction:
predicted sigma, signed residual. GC's equivalent is computed inline.
"""
import os
import argparse, csv, glob, json, math, re, sys
from pathlib import Path

import numpy as np
import pandas as pd

# ROOT is the eQuilibrator WORKING TREE, not this repository: these analyses read
# caches and fitted parameters that live there and are far too large to commit.
# It was __file__/../.. while this script lived in that tree; now it must be named
# explicitly. Overridable so the analysis is not pinned to one host.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))
# This repository, derived from the file location rather than hardcoded.
MSD = Path(__file__).resolve().parents[3]
DGP = Path(os.environ.get("DGPREDICTOR_DIR",
                          MSD.parent.parent / "dgpredictor_retrain"))
sys.path.insert(0, str(ROOT))
FOLDS, SEED, KCAL = 10, 20260901, 4.184
R_GAS = 8.31446261815324e-3


# ---------------------------------------------------------------- eQuilibrator
def calibrate_eq(cache, tecrdb, out):
    from component_contribution.trainer import Trainer
    from component_contribution.training_data import FullTrainingDataFactory
    from component_contribution.preprocessor import Preprocessor
    from equilibrator_cache import Q_
    from equilibrator_cache.reaction import Reaction
    from equilibrator_cache.api import create_compound_cache_from_sqlite_file

    ccache = create_compound_cache_from_sqlite_file(cache)
    td = FullTrainingDataFactory(ccache=ccache).make(tecrdb_path=tecrdb)
    S = td.stoichiometric_matrix                 # rows = Compound objects
    G = Trainer.group_incidence_matrix(td)
    compounds = list(S.index)
    b = np.array([float(x.m_as("kJ/mol")) if hasattr(x, "m_as") else float(x)
                  for x in td.standard_dg])
    w = np.array([float(x.m_as("dimensionless")) if hasattr(x, "m_as") else float(x)
                  for x in td.weight])
    ng = td.group_summary.shape[0]
    ncols = S.shape[1]
    print(f"eQ: {S.shape[0]} compounds x {ncols} training reactions", flush=True)

    Sid, Gid = S.copy(), G.copy()
    Sid.index = [c.id for c in compounds]
    Gid.index = [c.id for c in compounds]

    rng = np.random.default_rng(SEED)
    folds = np.array_split(rng.permutation(np.arange(ncols)), FOLDS)
    rec = []
    for i, fold in enumerate(folds, 1):
        train = np.setdiff1d(np.arange(ncols), fold)
        Ss = Sid.iloc[:, train]
        keep = (Ss != 0).any(axis=1)
        p = Trainer.train_from_matrices(
            Ss.loc[keep], Gid.loc[keep],
            pd.Series([Q_(float(x), "kJ/mol") for x in b[train]], index=Ss.columns),
            pd.Series(w[train], index=Ss.columns), ng)
        pre = Preprocessor(p)
        kept = set(np.asarray(Sid.index)[keep.values])
        n_ok = 0
        for col in fold:
            sparse = {compounds[r]: float(S.iloc[r, col])
                      for r in np.nonzero(S.iloc[:, col].values)[0]}
            if any(c.id not in kept for c in sparse):
                continue                      # not estimable from this fold
            try:
                mu, s_fin, s_inf, residual = pre.get_reaction_prediction(
                    Reaction(sparse))
            except Exception:
                continue
            if residual:
                continue                      # model declines; no sigma to score
            sigma = (np.linalg.norm(s_fin, 2)
                     + pre.RMSE_inf * np.linalg.norm(s_inf, 2))
            rec.append({"sigma_kcal": sigma / KCAL,
                        "residual_kcal": (mu - b[col]) / KCAL})
            n_ok += 1
        print(f"  fold {i}/{FOLDS}: scored {n_ok}/{len(fold)}", flush=True)
    _write(out, rec)
    return rec


# ---------------------------------------------------------------- dGPredictor
def calibrate_dg(out):
    from sklearn.linear_model import BayesianRidge
    z = np.load(DGP / "data/training_matrices.npz", allow_pickle=True)
    X, y = z["X_comb"], z["y"]
    print(f"dG: {X.shape[0]} training rows x {X.shape[1]} active fragments", flush=True)
    rng = np.random.default_rng(SEED)
    folds = np.array_split(rng.permutation(np.arange(X.shape[0])), FOLDS)
    rec = []
    for i, fold in enumerate(folds, 1):
        train = np.setdiff1d(np.arange(X.shape[0]), fold)
        m = BayesianRidge(tol=1e-6, fit_intercept=False, compute_score=False)
        m.fit(X[train], y[train])
        mu, sd = m.predict(X[fold], return_std=True)
        for p, s, t in zip(mu, sd, y[fold]):
            rec.append({"sigma_kcal": s / KCAL, "residual_kcal": (p - t) / KCAL})
        print(f"  fold {i}/{FOLDS}: scored {len(fold)}", flush=True)
    _write(out, rec)
    return rec


# ---------------------------------------------------------- Group contribution
def calibrate_gc(out):
    """GC never trains on TECRDB, so this is out-of-sample by construction.
    Restricted to reactions with no net bound-hydrogen change, where the
    Convention A / B difference vanishes and the comparison is meaningful."""
    k2c = {r["kegg_id"]: r["cpd_id"] for r in csv.DictReader(
        open(DGP / "data/training_compounds.tsv"), delimiter="\t")}
    gc, form = {}, {}
    for f in glob.glob(str(MSD / "Biochemistry/compound_*.json")):
        for e in json.load(open(f)):
            form[e["id"]] = e.get("formula")
            v = (e.get("thermodynamics") or {}).get("Group contribution")
            if v and abs(float(v[0])) < 1e7:
                gc[e["id"]] = (float(v[0]), float(v[1]))

    def hc(f):
        if not f:
            return None
        m = re.search(r"H(\d*)(?![a-z])", f)
        return 0 if not m else (int(m.group(1)) if m.group(1) else 1)

    rec = []
    for line in open(ROOT / "data/dgpredictor/TECRDB.tsv"):
        p = (line.rstrip("\n").split("\t") + [""] * 14)[:14]
        try:
            Kp = float(p[9]); T = float(p[10]) if p[10].strip() else 298.15
        except ValueError:
            continue
        if Kp <= 0 or not p[12].strip() or "=" not in p[6]:
            continue
        L, Rt = p[6].split("=", 1)
        st, ok = {}, True
        for side, sign in ((L, -1), (Rt, 1)):
            for term in re.split(r"\s+\+\s+", side.strip()):
                term = term.strip()
                if not term:
                    continue
                m = re.match(r"^(?:(\d+(?:\.\d+)?)\s+)?(C\d{5})$", term)
                if not m:
                    ok = False; break
                st[m.group(2)] = st.get(m.group(2), 0) + float(m.group(1) or 1) * sign
            if not ok:
                break
        if not ok:
            continue
        cpds = {}
        for k, n in st.items():
            c = k2c.get(k)
            if c is None or c not in gc:
                ok = False; break
            cpds[c] = cpds.get(c, 0) + n
        if not ok or not cpds:
            continue
        if any(hc(form.get(c)) is None for c in cpds):
            continue
        if sum(n * hc(form[c]) for c, n in cpds.items()) != 0:
            continue
        rec.append({"sigma_kcal": math.sqrt(sum((n * gc[c][1]) ** 2 for c, n in cpds.items())),
                    "residual_kcal": sum(n * gc[c][0] for c, n in cpds.items())
                                     + R_GAS * T * math.log(Kp) / KCAL})
    _write(out, rec)
    return rec


def _write(path, rec):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=["sigma_kcal", "residual_kcal"])
        w.writeheader()
        for r in rec:
            w.writerow({k: f"{v:.6g}" for k, v in r.items()})
    print(f"  wrote {path} ({len(rec)} rows)", flush=True)


def report(name, rec):
    if not rec:
        print(f"\n{name}: no rows"); return
    s = np.array([r["sigma_kcal"] for r in rec])
    e = np.array([abs(r["residual_kcal"]) for r in rec])
    fin = np.isfinite(s) & np.isfinite(e) & (s > 0)
    s, e = s[fin], e[fin]
    z = e / s
    print(f"\n{'='*78}\n{name}   n={len(s):,}")
    print(f"  median sigma {np.median(s):8.3f}   median |error| {np.median(e):8.3f} kcal/mol")
    print(f"  median |z| {np.median(z):6.3f}  (calibrated 0.674)   "
          f"frac|z|<1 {100*np.mean(z<1):5.1f}%  (calibrated 68%)")
    verdict = ("over-conservative -- error bars WIDER than the errors"
               if np.median(z) < 0.55 else
               "over-confident -- error bars NARROWER than the errors"
               if np.median(z) > 0.80 else "roughly calibrated")
    print(f"  -> {verdict}")
    print(f"\n  empirical scale for this source:")
    print(f"    {'reported sigma':<22}{'n':>7}{'median |error|':>16}{'p90 |error|':>13}")
    edges = [0, 0.5, 1, 2, 5, 10, 25, np.inf]
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (s >= lo) & (s < hi)
        if m.sum() < 5:
            continue
        lab = f"{lo:g} - {hi:g}" if np.isfinite(hi) else f"> {lo:g}"
        print(f"    {lab:<22}{int(m.sum()):>7}{np.median(e[m]):>16.3f}{np.percentile(e[m],90):>13.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=str(ROOT / "data/modelseed_cache_msd/compounds.sqlite"))
    ap.add_argument("--tecrdb", type=Path, default=ROOT / "data/opentecr/TECRDB_dedup.csv")
    ap.add_argument("--only", default=None, choices=("gc", "eq", "dg"))
    a = ap.parse_args()
    out = {}
    if a.only in (None, "gc"):
        out["Group contribution (vs measured, out-of-sample by construction)"] = \
            calibrate_gc(ROOT / "data/sigma_calibration_gc.tsv")
    if a.only in (None, "dg"):
        out["dGPredictor (10-fold held out)"] = \
            calibrate_dg(ROOT / "data/sigma_calibration_dg.tsv")
    if a.only in (None, "eq"):
        out["eQuilibrator (10-fold held out)"] = \
            calibrate_eq(a.cache, a.tecrdb, ROOT / "data/sigma_calibration_eq.tsv")
    for k, v in out.items():
        report(k, v)


if __name__ == "__main__":
    main()
