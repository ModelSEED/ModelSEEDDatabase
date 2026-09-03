#!/usr/bin/env python3
"""Grade every thermodynamic source, per reaction, gold / silver / bronze.

Sibling of ``optimize_thermo_source_assignment.py``, which picks ONE winning
source per reaction. This does not pick a winner -- it labels each source
independently, so on the same reaction eQuilibrator can be GOLD while Group
Contribution is BRONZE.

Four sources are graded:

    TECRDB                  experimental dG'^0. ALWAYS GOLD -- it is a
                            measurement, not a prediction. The one nuance is
                            the MATCH: the stereo_exact tier (full InChIKey)
                            is GOLD, the skeleton tier (connectivity only,
                            blind to anomers) is capped at SILVER because the
                            measurement may have been attached to the wrong
                            reaction. --tecrdb-skeleton-gold disables that.
    eQuilibrator            \
    dGPredictor-ModelSEED    > graded by the cascade below
    Group contribution      /

WHAT GRADES A PREDICTOR
-----------------------
1. ``p_ok = P(|dG_s - dG*| <= tau | sigma_s)``, tau = 2.0 kcal/mol -- the
   reversible band the direction cascade itself uses
   (reversibility_heuristics.py mmdeltag_band_heuristic). Fitted per source by
   isotonic regression on the same two-tier data as ehat (anchor = 802 TECRDB
   stereo-exact matches, weight 3; proxy = |source - trusted-sigma reference|,
   weight 1), but on the INDICATOR rather than the magnitude.

   Why not ehat: ehat is a good threshold and a weak ranking (rho 0.11-0.24 vs
   actual error), and on the Convention A Group Contribution rebuild it spans
   only 3.04-5.70 kcal/mol across the whole database -- a flat curve that
   cannot rank anything. A probability is also comparable across sources,
   which raw sigma is not (GC overstates its error 2.2x, dGPredictor-MS 1.5x,
   eQuilibrator understates 1.6x).

2. Direct measurement against TECRDB where it exists -- per source, so the
   three predictors are scored individually on the same reaction.

3. Cross-source behaviour, used ASYMMETRICALLY. Sources are fused by their
   calibrated error scale and scored by the Birge ratio R = sqrt(chi2/df) (the
   PDG scale factor for discrepant measurements), with per-source residual
   z_s = |dG_s - dG_fused| / ehat_s naming the outlier.

   Agreement between two fallible predictors is weak evidence -- they share
   lineage (eQuilibrator and Group Contribution both descend from group
   contribution) and 11% of concordant reactions are structural zeros where
   agreement is imposed by the stoichiometry. Disagreement is strong: someone
   is definitely wrong and z_s says who. Measured on TECRDB, letting
   corroboration promote to GOLD grew eQuilibrator's GOLD column 2,443 ->
   9,157 but diluted its measured guarantee 94% -> 90% within 2 kcal/mol.
   So: corroboration lifts BRONZE to SILVER and no further; being outvoted
   demotes one tier.

OUTPUTS (results/thermo_grades/)
    source_grades.tsv        long, one row per (reaction x source) -- primary
    source_grades_wide.tsv   one row per reaction, one grade column per source
    grade_calibration.json   fitted p_ok curves + validation tables
    grade_frontier.tsv       threshold sweep

CONSUMERS import ``load_grades()`` / ``recommended_energy_map()``.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
REPO_ROOT = Path(__file__).resolve().parents[3]
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", str(REPO_ROOT)))
MSDB_CODE = Path(os.environ.get("MSDB_CODE", str(REPO_ROOT)))
ANALYSIS_DIR = Path(os.environ.get("THERMO_GRADING_OUT",
                                   str(REPO_ROOT / "Biochemistry" / "Thermodynamics" / "SourceGrading")))
from optimize_thermo_source_assignment import (  # noqa: E402
    ANALYSIS_DIR, EQ_SENTINEL, SOURCES, fit_error_models, load_db, predict_error,
)

K = ["GC", "EQ", "DGPMS"]
LABEL = {"GC": "Group contribution", "EQ": "eQuilibrator",
         "DGPMS": "dGPredictor-ModelSEED", "TECRDB": "TECRDB"}
OUT = Path(os.environ.get("GRADES_OUT", str(ANALYSIS_DIR / "results" / "thermo_grades")))
TECRDB_CSV = Path(os.environ.get(
    "TECRDB_COMPARISON",
    str(REPO_ROOT / "Biochemistry" / "Thermodynamics" / "SourceGrading"
        / "tecrdb_vs_dgpredictor_modelseed.csv")))
RECONCILIATION = ANALYSIS_DIR / "results" / "eq_vs_dgpms" / "reconciliation.tsv"

TAU = 2.0            # kcal/mol -- the cascade's reversible half-band
P_GOLD = 0.90        # p_ok at or above this is self-certain
P_SILVER = 0.70
R_CORROB = 1.5       # Birge ratio at or below which sources corroborate
Z_CORROB = 1.0
R_OUTVOTE = 2.0      # Birge ratio above which the set is discrepant...
Z_OUTVOTE = 3.0      # ...and a source this far out is the one to blame
MEAS_GOLD = 1.0      # |dG_s - experiment| kcal/mol
MEAS_SILVER = 3.0
EHAT_FLOOR = 0.3     # kcal/mol; keeps 1/ehat^2 weights finite

GOLD, SILVER, BRONZE = 0, 1, 2
NAME = {GOLD: "GOLD", SILVER: "SILVER", BRONZE: "BRONZE"}
RNG = np.random.default_rng(20260812)


# --------------------------------------------------------------------- inputs
def load_tecrdb() -> pd.DataFrame:
    """Experimental dG'^0 per ModelSEED reaction, best match tier first."""
    t = pd.read_csv(TECRDB_CSV)
    t["tecrdb_dg"] = t["tecrdb_dG_kJ"] / 4.184
    t["tecrdb_sd"] = t["tecrdb_dG_sd_kJ"] / 4.184
    t = t[["modelseed_rxn", "match_tier", "tecrdb_dg", "tecrdb_sd",
           "n_measurements"]].rename(columns={"modelseed_rxn": "rxn"})
    # stereo_exact sorts before skeleton, so first() keeps the better tier
    return t.sort_values("match_tier").groupby("rxn", as_index=False).first()


def load_vetoes() -> set:
    """eQuilibrator reactions hit by the MetaNetX id collision bug.

    ``Retrieve_eQuilibrator_Reactions_Energies.py`` does ``lhs[mnx] = |coeff|``
    instead of accumulating, so two ModelSEED compounds sharing one MetaNetX id
    silently overwrite each other. Read from the reconciliation table if it has
    been built; absent, the veto is simply not applied (and said so).
    """
    if not RECONCILIATION.exists():
        print("  note: reconciliation.tsv absent, eQ MetaNetX-collision veto not applied")
        return set()
    d = pd.read_csv(RECONCILIATION, sep="\t", low_memory=False)
    if "eq_mnx_collision" not in d.columns:
        return set()
    return set(d.loc[d.eq_mnx_collision == 1, "rxn"])


# ---------------------------------------------------------------- calibration
def fit_p_ok(train: pd.DataFrame, db: pd.DataFrame, tau: float = TAU) -> dict:
    """sigma -> P(|error| <= tau), per source, isotonic and DECREASING.

    Same two-tier anchor/proxy construction as ``fit_error_models``; see that
    docstring for why gold data alone cannot calibrate the sigma range the
    model has to work over.
    """
    from sklearn.isotonic import IsotonicRegression
    from optimize_thermo_source_assignment import (
        EQ_TRUSTED_SIGMA, DGPMS_TRUSTED_SIGMA)
    proxy = {"GC": "EQ", "DGPMS": "EQ", "EQ": "DGPMS"}
    trusted = {"EQ": EQ_TRUSTED_SIGMA, "DGPMS": DGPMS_TRUSTED_SIGMA}
    models = {}
    for k in K:
        m = train[f"dg_{k}"].notna() & train[f"sig_{k}"].notna()
        if k == "EQ":
            m &= train[f"sig_{k}"] <= EQ_SENTINEL
        sub = train[m]
        xs = list(sub[f"sig_{k}"].to_numpy(float))
        ys = list(((sub[f"dg_{k}"] - sub["tecrdb_dg"]).abs() <= tau).astype(float))
        w = [3.0] * len(xs)
        n_gold = len(xs)

        ref = proxy[k]
        ok = (db[f"dg_{k}"].notna() & db[f"sig_{k}"].notna()
              & db[f"dg_{ref}"].notna() & (db[f"sig_{ref}"] <= trusted[ref]))
        if k == "EQ":
            ok &= db[f"sig_{k}"] <= EQ_SENTINEL
        if k == "DGPMS":
            ok &= db["is_quinone"] == 0
        sv = db[ok]
        n_silver = len(sv)
        xs += list(sv[f"sig_{k}"].to_numpy(float))
        ys += list(((sv[f"dg_{k}"] - sv[f"dg_{ref}"]).abs() <= tau).astype(float))
        w += [1.0] * n_silver

        iso = IsotonicRegression(increasing=False, out_of_bounds="clip").fit(
            np.array(xs), np.array(ys), sample_weight=np.array(w))
        models[k] = {"x": iso.X_thresholds_.tolist(), "y": iso.y_thresholds_.tolist(),
                     "n_anchor": n_gold, "n_proxy": n_silver,
                     "anchor_frac_within_tau": float(np.mean(ys[:n_gold])) if n_gold else np.nan}
    return models


def predict_p_ok(db: pd.DataFrame, models: dict, veto_eq: set) -> pd.DataFrame:
    """p_ok per reaction per source; NaN where the source is absent or vetoed."""
    out = pd.DataFrame(index=db.index)
    for k in K:
        mdl = models[k]
        p = pd.Series(np.nan, index=db.index)
        avail = db[f"dg_{k}"].notna() & db[f"sig_{k}"].notna()
        p[avail] = np.interp(db.loc[avail, f"sig_{k}"].to_numpy(float),
                             mdl["x"], mdl["y"], left=mdl["y"][0], right=mdl["y"][-1])
        if k == "EQ":
            p[db["sig_EQ"] > EQ_SENTINEL] = np.nan       # source disclaims it
            p[db["rxn"].isin(veto_eq)] = np.nan          # MetaNetX collision
        if k == "DGPMS":
            p[db["is_quinone"] == 1] = np.nan            # 52.8% sign-wrong
        out[f"p_{k}"] = p
    return out


# ------------------------------------------------------------------- fusion
def fuse(db: pd.DataFrame, ehat: pd.DataFrame, p_ok: pd.DataFrame) -> pd.DataFrame:
    """Precision-weighted combination + Birge ratio + per-source residual.

    dG_fused is an internal construct -- a reference point for computing z_s --
    not a recommended value. Weights use the CALIBRATED error scale ehat, not
    raw sigma, because the three sigma scales are not commensurable.
    """
    dg = np.array([db[f"dg_{k}"].to_numpy(float) for k in K])
    tau = np.array([np.maximum(ehat[f"ehat_{k}"].to_numpy(float), EHAT_FLOOR) for k in K])
    usable = np.array([p_ok[f"p_{k}"].notna().to_numpy() for k in K])
    ok = np.isfinite(dg) & np.isfinite(tau) & usable

    w = np.where(ok, 1.0 / np.where(ok, tau, 1.0) ** 2, 0.0)
    wsum = w.sum(0)
    n = ok.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        fused = np.where(wsum > 0, (np.where(ok, dg, 0.0) * w).sum(0) / np.where(wsum > 0, wsum, np.nan), np.nan)
        chi2 = (np.where(ok, (dg - fused) ** 2, 0.0) * w).sum(0)
        R = np.where(n >= 2, np.sqrt(chi2 / np.maximum(n - 1, 1)), np.nan)
        tau_fused = np.sqrt(1.0 / np.where(wsum > 0, wsum, np.nan)) * np.maximum(np.nan_to_num(R, nan=1.0), 1.0)

    out = pd.DataFrame({"n_src": n, "dg_fused": fused, "birge": R,
                        "tau_fused": tau_fused}, index=db.index)
    # structural zero: every source says ~0, or the reaction is transport.
    # Their agreement is imposed by the chemistry, not earned.
    near_zero = np.where(ok, np.abs(dg) < 0.5, True).all(0) & (n >= 2)
    out["struct_zero"] = near_zero | (db["is_transport"] == 1).to_numpy()
    for i, k in enumerate(K):
        out[f"z_{k}"] = np.where(ok[i], np.abs(dg[i] - fused) / tau[i], np.nan)
    return out


# -------------------------------------------------------------------- grading
def grade_predictors(db, p_ok, fus, tec, disable_measured=False):
    """The cascade, applied independently per source. Returns (grades, reasons)."""
    grades, reasons = {}, {}
    meas_ok = tec["match_tier"].eq("stereo_exact") if "match_tier" in tec else pd.Series(False, index=db.index)
    corrob = (fus.n_src >= 2) & (~fus.struct_zero) & (fus.birge <= R_CORROB)
    outvote = (fus.n_src >= 2) & (fus.birge > R_OUTVOTE)

    for k in K:
        p = p_ok[f"p_{k}"]
        z = fus[f"z_{k}"]
        g = pd.Series(BRONZE, index=db.index)
        r = pd.Series("uncorroborated", index=db.index)

        # 2. base grade from the source's own calibrated confidence
        sel = p >= P_SILVER
        g[sel], r[sel] = SILVER, "self-confident"
        sel = p >= P_GOLD
        g[sel], r[sel] = GOLD, "self-certain"

        # 3. corroboration is a FLOOR, never a promotion to gold
        sel = corrob & (z <= Z_CORROB) & (g == BRONZE)
        g[sel], r[sel] = SILVER, "corroborated"

        # 4. being the outlier in a discrepant set costs one tier
        sel = outvote & (z > Z_OUTVOTE)
        g[sel] = np.minimum(g[sel] + 1, BRONZE)
        r[sel] = "outvoted"

        # 1. measurement overrides everything (applied last so it wins)
        if not disable_measured:
            err = (db[f"dg_{k}"] - tec["tecrdb_dg"]).abs()
            sel = meas_ok & err.notna()
            g[sel & (err <= MEAS_GOLD)] = GOLD
            g[sel & (err > MEAS_GOLD) & (err <= MEAS_SILVER)] = SILVER
            g[sel & (err > MEAS_SILVER)] = BRONZE
            r[sel] = "measured"

        g[p.isna()] = np.nan
        r[p.isna()] = "ungraded"
        grades[k], reasons[k] = g, r
    return grades, reasons


def grade_tecrdb(tec: pd.DataFrame, skeleton_gold: bool) -> tuple:
    g = pd.Series(np.nan, index=tec.index)
    r = pd.Series("ungraded", index=tec.index)
    have = tec["tecrdb_dg"].notna()
    g[have], r[have] = GOLD, "measured"
    if not skeleton_gold:
        skel = have & tec["match_tier"].eq("skeleton")
        g[skel], r[skel] = SILVER, "measured (skeleton match)"
    return g, r


# ----------------------------------------------------------------- validation
def validate(db, p_ok, fus, tec) -> list:
    """Regrade with the measurement override OFF, then score against it.

    IN-SAMPLE with respect to the calibration: the curves were fitted on all 802
    anchors and this scores those same 802. Disabling the measurement override
    removes the *direct* use of the label, not the indirect one. Use
    ``validate_cv`` for the leak-free number; the two agree to 0.04 kcal/mol
    (§ the docs), because the anchor is only 19-37% of the fitting weight and
    isotonic pooling makes one point out of 802 nearly invisible.
    """
    grades, _ = grade_predictors(db, p_ok, fus, tec, disable_measured=True)
    m0 = tec["match_tier"].eq("stereo_exact")
    rows = []
    for k in K:
        err = (db[f"dg_{k}"] - tec["tecrdb_dg"]).abs()
        for v in (GOLD, SILVER, BRONZE):
            m = m0 & (grades[k] == v) & err.notna()
            if not m.sum():
                rows.append({"source": LABEL[k], "grade": NAME[v], "n": 0})
                continue
            e = err[m]
            rows.append({"source": LABEL[k], "grade": NAME[v], "n": int(m.sum()),
                         "median_abs_err": float(e.median()),
                         "mean_abs_err": float(e.mean()),
                         "frac_within_1": float((e <= 1).mean()),
                         "frac_within_2": float((e <= 2).mean()),
                         "p90_abs_err": float(e.quantile(0.9))})
    return rows


def validate_cv(db, tec, veto_eq, n_folds: int = 5, n_reps: int = 4,
                seed: int = 11) -> list:
    """Leak-free validation: refit the calibration on each training fold.

    ``validate`` scores the same 802 anchors the curves were fitted on. Here the
    anchors are split into folds; for each fold the magnitude and probability
    curves are refitted on the *other* folds only, test reactions are also kept
    out of the proxy tier (matching the guard in
    ``optimize_thermo_source_assignment.main``), and only then is the held-out
    fold graded and scored.

    Slower -- it refits everything ``n_folds * n_reps`` times -- so it is not on
    the default path. ``--cv`` runs it.
    """
    m0 = tec["match_tier"].eq("stereo_exact")
    anchor_idx = db.index[m0.to_numpy() & db["tecrdb_dg"].notna()]
    rng = np.random.default_rng(seed)
    pooled = {(k, v): [] for k in K for v in (GOLD, SILVER, BRONZE)}

    for _ in range(n_reps):
        folds = np.array_split(rng.permutation(anchor_idx), n_folds)
        for f in range(n_folds):
            te = folds[f]
            tr = np.concatenate([folds[j] for j in range(n_folds) if j != f])
            db_tr = db[~db.rxn.isin(set(db.loc[te, "rxn"]))]   # test out of proxy too
            eh = predict_error(db, fit_error_models(
                db.loc[tr].rename(columns={"tecrdb_dg": "tecrdb"}), db_tr))
            pk = predict_p_ok(db, fit_p_ok(db.loc[tr], db_tr), veto_eq)
            gr, _ = grade_predictors(db, pk, fuse(db, eh, pk), tec, disable_measured=True)
            for k in K:
                err = (db[f"dg_{k}"] - db["tecrdb_dg"]).abs()
                for v in (GOLD, SILVER, BRONZE):
                    m = pd.Series(False, index=db.index)
                    m.loc[te] = True
                    m &= (gr[k] == v) & err.notna()
                    pooled[(k, v)].extend(err[m].tolist())

    rows = []
    for k in K:
        for v in (GOLD, SILVER, BRONZE):
            e = np.array(pooled[(k, v)])
            if not len(e):
                rows.append({"source": LABEL[k], "grade": NAME[v], "n": 0})
                continue
            rows.append({"source": LABEL[k], "grade": NAME[v], "n": int(len(e)),
                         "median_abs_err": float(np.median(e)),
                         "mean_abs_err": float(e.mean()),
                         "frac_within_1": float((e <= 1).mean()),
                         "frac_within_2": float((e <= 2).mean()),
                         "p90_abs_err": float(np.percentile(e, 90)),
                         "note": f"{n_folds}-fold x {n_reps} reps; n counts each "
                                 f"anchor once per rep"})
    return rows


# --------------------------------------------------------------------- output
def build(skeleton_gold: bool = False, heldout: bool = False) -> tuple:
    """``heldout=True`` grades the three predictors with the TECRDB measurement
    override DISABLED and omits TECRDB as a source.

    This exists so the graded map can be scored against TECRDB without
    circularity: the ordinary grades use the measurement, so a graded direction
    map that includes TECRDB necessarily reproduces TECRDB perfectly. The
    held-out grades never see it.
    """
    db = load_db().reset_index(drop=True)
    tec_raw = load_tecrdb()
    tec = db[["rxn"]].merge(tec_raw, on="rxn", how="left")
    tec.index = db.index
    db["tecrdb_dg"] = tec["tecrdb_dg"]

    anchor = db[tec["match_tier"].eq("stereo_exact").to_numpy() & db["tecrdb_dg"].notna()]
    print(f"non-EMPTY reactions {len(db)};  TECRDB stereo-exact anchor {len(anchor)}, "
          f"skeleton {int(tec['match_tier'].eq('skeleton').sum())}")

    veto_eq = load_vetoes()
    ehat = predict_error(db, fit_error_models(anchor.rename(columns={"tecrdb_dg": "tecrdb"}), db))
    pmods = fit_p_ok(anchor, db)
    p_ok = predict_p_ok(db, pmods, veto_eq)
    fus = fuse(db, ehat, p_ok)

    print("\ncalibration  p_ok = P(|error| <= %.1f kcal/mol | sigma):" % TAU)
    for k in K:
        m = pmods[k]
        print(f"  {LABEL[k]:22s} anchor n={m['n_anchor']:4d} proxy n={m['n_proxy']:6d}  "
              f"p_ok {m['y'][0]:.3f} -> {m['y'][-1]:.3f}  "
              f"anchor frac within tau {m['anchor_frac_within_tau']:.3f}")

    val = validate(db, p_ok, fus, tec)
    print("\nvalidation on TECRDB stereo-exact, measurement override DISABLED:")
    for r in val:
        if r["n"]:
            print(f"  {r['source']:22s} {r['grade']:6s} n={r['n']:4d}  "
                  f"median={r['median_abs_err']:6.2f}  <=2: {r['frac_within_2']:4.0%}")
        else:
            print(f"  {r['source']:22s} {r['grade']:6s} n=   0")

    grades, reasons = grade_predictors(db, p_ok, fus, tec, disable_measured=heldout)
    tg, tr = grade_tecrdb(tec, skeleton_gold)

    long_rows = []
    for k in K:
        long_rows.append(pd.DataFrame({
            "rxn": db.rxn, "name": db.name, "ec": db.ec, "source": LABEL[k],
            "dg": db[f"dg_{k}"], "sigma": db[f"sig_{k}"], "operator": db.get(f"op_{k}"),
            "ehat": ehat[f"ehat_{k}"], "p_ok": p_ok[f"p_{k}"], "z": fus[f"z_{k}"],
            "birge": fus.birge, "n_src": fus.n_src, "struct_zero": fus.struct_zero,
            "grade": grades[k].map(NAME), "reason": reasons[k],
            "n_anchor": pmods[k]["n_anchor"], "n_proxy": pmods[k]["n_proxy"]}))
    if not heldout:
        long_rows.append(pd.DataFrame({
            "rxn": db.rxn, "name": db.name, "ec": db.ec, "source": "TECRDB",
            "dg": tec.tecrdb_dg, "sigma": tec.tecrdb_sd, "operator": None,
            "ehat": np.nan, "p_ok": np.nan, "z": np.nan, "birge": fus.birge,
            "n_src": fus.n_src, "struct_zero": fus.struct_zero,
            "grade": tg.map(NAME), "reason": tr,
            "n_anchor": np.nan, "n_proxy": np.nan}))
    long = pd.concat(long_rows, ignore_index=True)
    long = long[long.grade.notna()].copy()

    wide = db[["rxn", "name", "ec", "status", "is_transport"]].copy()
    wide["n_src"] = fus.n_src
    wide["birge"] = fus.birge
    wide["struct_zero"] = fus.struct_zero
    for k in K:
        wide[f"dg_{k}"] = db[f"dg_{k}"]
        wide[f"grade_{k}"] = grades[k].map(NAME)
    wide["dg_TECRDB"] = tec.tecrdb_dg
    wide["grade_TECRDB"] = tg.map(NAME)
    ranks = pd.concat([tg.rename("TECRDB")] + [grades[k].rename(k) for k in K], axis=1)
    wide["best_grade"] = ranks.min(axis=1).map(NAME)
    wide["best_source"] = ranks.idxmin(axis=1).where(ranks.notna().any(axis=1))
    calib = {"tau": TAU, "thresholds": {"p_gold": P_GOLD, "p_silver": P_SILVER,
                                        "r_corrob": R_CORROB, "z_corrob": Z_CORROB,
                                        "r_outvote": R_OUTVOTE, "z_outvote": Z_OUTVOTE,
                                        "meas_gold": MEAS_GOLD, "meas_silver": MEAS_SILVER},
             "tecrdb_skeleton_gold": skeleton_gold,
             "msdb_root": str(MSDB_ROOT),
             "p_ok_models": pmods, "validation": val,
             "vetoes": {"eq_sentinel": int((db.sig_EQ > EQ_SENTINEL).sum()),
                        "eq_mnx_collision": len(veto_eq),
                        "dgpms_quinone": int(((db.is_quinone == 1) & db.dg_DGPMS.notna()).sum())}}
    calib["validation_note"] = (
        "validation[] is IN-SAMPLE w.r.t. the calibration: the curves were fitted "
        "on all 802 anchors and this scores those same 802, with only the "
        "measurement override disabled. validation_cv[], when present, refits per "
        "fold and is leak-free.")
    return long, wide, calib, (db, tec, veto_eq)


def frontier(long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for src, g in long.groupby("source"):
        c = g.grade.value_counts()
        rows.append({"source": src, "n_graded": len(g),
                     **{f"n_{n.lower()}": int(c.get(n, 0)) for n in ("GOLD", "SILVER", "BRONZE")},
                     **{f"n_reason_{r}": int(v) for r, v in g.reason.value_counts().items()}})
    return pd.DataFrame(rows).fillna(0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tecrdb-skeleton-gold", action="store_true",
                    help="grade skeleton-tier TECRDB matches GOLD instead of SILVER")
    ap.add_argument("--cv", action="store_true",
                    help="also run the leak-free cross-validated validation "
                         "(refits the calibration per fold; several minutes)")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    long, wide, calib, ctx = build(skeleton_gold=args.tecrdb_skeleton_gold)

    print("\n=== grades ===")
    tab = pd.crosstab(long.source, long.grade).reindex(
        index=[LABEL["TECRDB"], LABEL["EQ"], LABEL["DGPMS"], LABEL["GC"]],
        columns=["GOLD", "SILVER", "BRONZE"]).fillna(0).astype(int)
    print(tab.to_string())
    print("\nper reaction:")
    bg = wide.best_grade.value_counts()
    print(f"  any source {int(wide.best_grade.notna().sum())}   "
          + "   ".join(f"best={n}: {int(bg.get(n, 0))}" for n in ("GOLD", "SILVER", "BRONZE")))

    if args.cv:
        db_, tec_, veto_ = ctx
        print("\nleak-free cross-validated validation (refits per fold)...")
        calib["validation_cv"] = validate_cv(db_, tec_, veto_)
        for r in calib["validation_cv"]:
            if r["n"]:
                print(f"  {r['source']:22s} {r['grade']:6s} n={r['n']:5d}  "
                      f"median={r['median_abs_err']:6.2f}  <=2: {r['frac_within_2']:4.0%}")

    long.to_csv(OUT / "source_grades.tsv", sep="\t", index=False, float_format="%.4f")
    wide.to_csv(OUT / "source_grades_wide.tsv", sep="\t", index=False, float_format="%.4f")
    frontier(long).to_csv(OUT / "grade_frontier.tsv", sep="\t", index=False)
    (OUT / "grade_calibration.json").write_text(json.dumps(calib, indent=1))
    print(f"\nwrote {OUT}/source_grades.tsv ({len(long)} rows), source_grades_wide.tsv, "
          "grade_frontier.tsv, grade_calibration.json")

    # held-out grades: no TECRDB source, no measurement override. Needed to
    # score a graded direction map against TECRDB without circularity.
    ho, _, _, _ = build(skeleton_gold=args.tecrdb_skeleton_gold, heldout=True)
    ho.to_csv(OUT / "source_grades_heldout.tsv", sep="\t", index=False, float_format="%.4f")
    print(f"wrote {OUT}/source_grades_heldout.tsv ({len(ho)} rows)")


# ------------------------------------------------------------------ consumers
def load_grades(wide: bool = False) -> pd.DataFrame:
    path = OUT / ("source_grades_wide.tsv" if wide else "source_grades.tsv")
    if not path.exists():
        raise FileNotFoundError(f"{path} missing; run grade_thermo_sources.py")
    return pd.read_csv(path, sep="\t", low_memory=False)


GRADE_ORDER = {"GOLD": 0, "SILVER": 1, "BRONZE": 2}


def recommended_energy_map(min_grade: str = "BRONZE", heldout: bool = False) -> dict:
    """``{rxn_id: (dg, dge, source_label)}`` -- the best-graded source per reaction.

    Ranked by grade, then by p_ok descending (TECRDB, having no p_ok, always
    sorts first within GOLD because it is scored as a measurement). Reactions
    whose best grade is worse than ``min_grade`` are omitted entirely, so a
    downstream consumer sees them as "no data" rather than as a bad number.
    """
    floor = GRADE_ORDER[min_grade]
    path = OUT / ("source_grades_heldout.tsv" if heldout else "source_grades.tsv")
    if not path.exists():
        raise FileNotFoundError(f"{path} missing; run grade_thermo_sources.py")
    g = pd.read_csv(path, sep="\t", low_memory=False)
    g = g[g.grade.map(GRADE_ORDER) <= floor].copy()
    g["_g"] = g.grade.map(GRADE_ORDER)
    g["_p"] = np.where(g.source == "TECRDB", 2.0, g.p_ok.fillna(-1.0))
    g = g.sort_values(["rxn", "_g", "_p"], ascending=[True, True, False])
    best = g.groupby("rxn", as_index=False).first()
    out = {}
    for r in best.itertuples():
        if not np.isfinite(r.dg):
            continue
        sd = r.sigma if np.isfinite(r.sigma) else 0.0
        out[r.rxn] = (float(r.dg), float(sd), r.source)
    return out


if __name__ == "__main__":
    main()
