#!/usr/bin/env python3
"""Choose WHICH thermodynamic source to use for each ModelSEED reaction, so that
coverage is maximised at a stated error tolerance.

This is a different problem from optimize_thermo_consensus.py, and a bigger one.
That script picks a SUBSET of reactions where eQuilibrator and dGPredictor
agree -- a consensus set, ceiling ~3k reactions, because it needs both sources
present and concordant. This script makes a PER-REACTION SOURCE ASSIGNMENT, so a
reaction with only one usable source still contributes. Ceiling: 32,466
reactions (58% of the database), roughly 10x.

    for each reaction i:   s*(i) = argmin_s  ehat(i, s)   over available sources
                           keep i  iff  min_s ehat(i, s) <= E*
    maximise               coverage = |{i kept}|

Given ehat this is solved exactly by inspection -- reactions do not interact, so
there is no combinatorial search. ALL of the difficulty is in ehat, the expected
absolute error of source s on reaction i.

WHY ehat NEEDS GROUND TRUTH
---------------------------
Cross-source disagreement gives JOINT error, not per-source error: if two
sources differ by 20 kcal/mol you know one is wrong, not which. Splitting them
requires an external reference. We use TECRDB (NIST experimental dG'^0), matched
to ModelSEED reactions by the SMILES->InChIKey multiset pipeline in
/scratch/ctaylor/dgpredictor_tecrdb (1,550 matches, of which 802 are
``stereo_exact`` -- the tier that distinguishes anomers and D/L pairs, and the
only tier used for fitting here).

On that reference, median |error| is eQuilibrator 0.45, dGPredictor-ModelSEED
0.47, Group Contribution 1.60 kcal/mol -- so no source dominates and the
assignment is worth making.

CALIBRATION
-----------
For each source we regress observed |error vs TECRDB| on that source's OWN
reported sigma, using isotonic regression: monotone (more self-reported
uncertainty must not predict less error) and non-parametric (no functional form
imposed on a relationship we have no theory for). ehat(i,s) is then that
source's fitted curve evaluated at its sigma for reaction i, which extends to
every reaction the source covers, not just the matched ones.

Two hard overrides sit on top, both established elsewhere in this analysis and
neither visible to a sigma-only model:
  * eQuilibrator sentinel uncertainties (>100 kcal/mol) -- the source declaring
    it has no estimate. Never assigned.
  * dGPredictor on the quinone / quinol couple -- 52.8% sign-wrong, and
    self-flagged with a median sigma of 80 (report section 2). Never assigned.

VALIDATION
----------
Fitted on a TECRDB train split, scored on held-out TECRDB reactions against four
baselines including the incumbent -- dev's
Promote_Reaction_Thermodynamics_to_Canonical.py priority, eQuilibrator then
Group Contribution then the ML tier, lowest reported error within a tier.

OUTPUTS (results/eq_vs_dgpms/)
    source_assignment.tsv        per reaction: chosen source, merged dG, ehat
    source_assignment_frontier.tsv   coverage vs error tolerance
    source_assignment_models.json    fitted calibration + validation

CONSUMERS import ``load_assignment()``.
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --- paths -------------------------------------------------------------
# Relocated into ModelSEEDDatabase 2026-09-03. Previously this lived in
# core_models_analysis and reached back into the database through MSDB_ROOT and
# MSDB_CODE; both now resolve to this repository, so the cross-repo hop is gone
# and `reversibility_heuristics` is a sibling import rather than a sys.path
# injection against another checkout.
#
# Only the grading and recommendation stages moved. The core-model FBA stages
# stayed behind: they need 5,683 models and cobra/GLPK, which is analysis rather
# than database.
REPO_ROOT = Path(__file__).resolve().parents[3]
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", str(REPO_ROOT)))
MSDB_CODE = Path(os.environ.get("MSDB_CODE", str(REPO_ROOT)))
BIOCHEM = MSDB_ROOT / "Biochemistry"
ANALYSIS_DIR = Path(os.environ.get("THERMO_GRADING_OUT",
                                   str(REPO_ROOT / "Biochemistry" / "Thermodynamics"
                                       / "SourceGrading")))

sys.path.insert(0, str(Path(__file__).resolve().parent))
from organic_reaction_types import QUINONE_RE  # noqa: E402

# NOTE: the earlier snapshot /scratch/ctaylor/tmp/devsnap (dev @ 34992d39) was
# deleted 2026-08-12. devsnap2 is dev @ 49563c6f: eQuilibrator and
# dGPredictor-ModelSEED are byte-identical to it, Group Contribution is the
# Convention A rebuild (ad34d6ab) -- 53% of values changed, coverage +1,501.
# Re-running with these defaults therefore refits GC against Convention A and
# would overwrite results/eq_vs_dgpms/, which was fitted on the OLD GC and is
# what EQUILIBRATOR_VS_DGPREDICTOR_MODELSEED.md quotes. The Convention A refit
# already exists as results/eq_vs_dgpms_gcA/ -- set EQDGP_OUT to keep them apart.
MSDB_ROOT = MSDB_ROOT
ANALYSIS_DIR = ANALYSIS_DIR
OUT = Path(os.environ.get("EQDGP_OUT", str(ANALYSIS_DIR / "results" / "eq_vs_dgpms")))
TECRDB = Path(os.environ.get(
    "TECRDB_COMPARISON",
    str(REPO_ROOT / "Biochemistry" / "Thermodynamics" / "SourceGrading"
        / "tecrdb_vs_dgpredictor_modelseed.csv")))
BIOCHEM = MSDB_ROOT / "Biochemistry"

SOURCES = {"Group contribution": "GC", "eQuilibrator": "EQ",
           "dGPredictor-ModelSEED": "DGPMS"}
EQ_SENTINEL = 100.0        # kcal/mol; eQuilibrator's "no estimate" marker
TOLERANCE = 2.0            # kcal/mol expected error, the shipped operating point
RNG = np.random.default_rng(20260806)

MODELS_JSON = OUT / "source_assignment_models.json"
ASSIGN_TSV = OUT / "source_assignment.tsv"


# ------------------------------------------------------------------- loading
def load_db() -> pd.DataFrame:
    """One row per non-EMPTY reaction: each source's dG and sigma, quinone flag."""
    cpd_name = {}
    for path in sorted(glob.glob(str(BIOCHEM / "compound_*.json"))):
        for e in json.load(open(path)):
            cpd_name[e["id"]] = str(e.get("name") or "")

    rows = []
    for path in sorted(glob.glob(str(BIOCHEM / "reaction_*.json"))):
        for e in json.load(open(path)):
            if e.get("status") == "EMPTY":
                continue
            th = e.get("thermodynamics") or {}
            r = {"rxn": e["id"], "name": e.get("name", ""),
                 "ec": ";".join(e.get("ec_numbers") or []),
                 "status": e.get("status", ""),
                 "is_transport": e.get("is_transport", 0)}
            for lbl, key in SOURCES.items():
                t = th.get(lbl)
                if t and len(t) > 2 and t[2] not in (None, "?"):
                    try:
                        r[f"dg_{key}"] = float(t[0])
                        r[f"sig_{key}"] = abs(float(t[1]))
                        r[f"op_{key}"] = t[2]
                    except (TypeError, ValueError):
                        pass
            names = [cpd_name.get(s["compound"], "")
                     for s in (e.get("stoichiometry") or [])]
            r["is_quinone"] = int(any(QUINONE_RE.search(n) for n in names))
            rows.append(r)
    return pd.DataFrame(rows)


def load_truth(db: pd.DataFrame) -> pd.DataFrame:
    """TECRDB experimental dG'^0, stereo_exact tier only, joined to the snapshot."""
    t = pd.read_csv(TECRDB)
    t = t[t.match_tier == "stereo_exact"].copy()
    t["tecrdb"] = t["tecrdb_dG_kJ"] / 4.184
    t = t[["modelseed_rxn", "tecrdb", "n_measurements"]].rename(
        columns={"modelseed_rxn": "rxn"})
    t = t.groupby("rxn", as_index=False).first()
    return db.merge(t, on="rxn", how="inner")


# --------------------------------------------------------------- calibration
# eQuilibrator sigma below which TECRDB shows it is accurate to ~0.45 kcal/mol.
# This is its TECRDB p90, so it is the range where the gold data actually
# constrains it -- see PROXY REFERENCE below.
EQ_TRUSTED_SIGMA = 0.70
DGPMS_TRUSTED_SIGMA = 1.22


def fit_error_models(train: pd.DataFrame, db: pd.DataFrame | None = None) -> dict:
    """sigma -> expected |error| per source, by isotonic regression on two tiers.

    Monotone (a source reporting more uncertainty must not be predicted more
    accurate) and non-parametric (no theory for the shape).

    WHY TWO TIERS. TECRDB only covers well-measured central metabolism, which is
    exactly the LOW-sigma regime, so gold data alone cannot calibrate the range
    the model must actually work over:

        dGPredictor-MS sigma   TECRDB p50 0.91, p90 1.22, max 21.6
                               database p50 21.17, p90 52.89, max 2039

    75.6% of database reactions for dGPredictor (43.4% eQuilibrator, 27.8% Group
    Contribution) lie beyond the TECRDB p90. Fitting on gold alone and clipping
    would assign them the error learned at sigma ~ 1.2 -- underestimating error
    precisely where the source is least reliable, which is the opposite of what
    a safety filter must do.

    PROXY REFERENCE. TECRDB establishes that eQuilibrator with sigma <= 0.70 is
    accurate to a median 0.45 kcal/mol. That earns it the right to stand in as a
    reference where gold data runs out: for the other sources the silver target
    is |source - eQuilibrator| on reactions where eQuilibrator is in its trusted
    band. eQuilibrator itself is scored the same way against low-sigma
    dGPredictor. The silver target is a STAND-IN, not a measurement. Measured
    against the anchor it is approximately unbiased -- median difference +/-0.01
    kcal/mol, and it exceeds the real error only ~half the time -- but noisier
    (Spearman 0.43-0.84). So gold is weighted 3x for lower variance, not to
    correct a bias, and the tier each reaction was fitted from is recorded.
    """
    from sklearn.isotonic import IsotonicRegression
    models = {}
    proxy = {"GC": "EQ", "DGPMS": "EQ", "EQ": "DGPMS"}
    trusted = {"EQ": EQ_TRUSTED_SIGMA, "DGPMS": DGPMS_TRUSTED_SIGMA}
    for k in SOURCES.values():
        m = train[f"dg_{k}"].notna() & train[f"sig_{k}"].notna()
        if k == "EQ":
            m &= train[f"sig_{k}"] <= EQ_SENTINEL
        sub = train[m]
        xs = list(sub[f"sig_{k}"].to_numpy(float))
        ys = list((sub[f"dg_{k}"] - sub["tecrdb"]).abs().to_numpy(float))
        w = [3.0] * len(xs)
        n_gold = len(xs)

        n_silver = 0
        ref = proxy[k]
        if db is not None and f"dg_{ref}" in db:
            ok = (db[f"dg_{k}"].notna() & db[f"sig_{k}"].notna()
                  & db[f"dg_{ref}"].notna() & (db[f"sig_{ref}"] <= trusted[ref]))
            if k == "EQ":
                ok &= db[f"sig_{k}"] <= EQ_SENTINEL
            if k == "DGPMS":
                ok &= db["is_quinone"] == 0      # excluded by override anyway
            sv = db[ok]
            n_silver = len(sv)
            xs += list(sv[f"sig_{k}"].to_numpy(float))
            ys += list((sv[f"dg_{k}"] - sv[f"dg_{ref}"]).abs().to_numpy(float))
            w += [1.0] * n_silver

        if n_gold + n_silver < 30:
            models[k] = {"kind": "const", "n_gold": n_gold, "n_silver": n_silver,
                         "value": float(np.median(ys)) if ys else 10.0}
            continue
        xs, ys, w = np.array(xs), np.array(ys), np.array(w)
        iso = IsotonicRegression(increasing=True, out_of_bounds="clip").fit(xs, ys,
                                                                           sample_weight=w)
        gold_x = sub[f"sig_{k}"].to_numpy(float)
        models[k] = {"kind": "isotonic", "n_gold": n_gold, "n_silver": n_silver,
                     "x": iso.X_thresholds_.tolist(), "y": iso.y_thresholds_.tolist(),
                     "gold_median_err": float(np.median(
                         (sub[f"dg_{k}"] - sub["tecrdb"]).abs())) if n_gold else np.nan,
                     "gold_sigma_p90": float(np.percentile(gold_x, 90)) if n_gold else np.nan,
                     "spearman_sigma_vs_err": float(
                         pd.Series(xs).corr(pd.Series(ys), method="spearman"))}
    return models


def predict_error(db: pd.DataFrame, models: dict) -> pd.DataFrame:
    """ehat(i, s) for every reaction and every source it covers."""
    out = pd.DataFrame(index=db.index)
    for k, mdl in models.items():
        sig = db.get(f"sig_{k}")
        e = pd.Series(np.nan, index=db.index)
        avail = db[f"dg_{k}"].notna() if f"dg_{k}" in db else pd.Series(False, index=db.index)
        if mdl["kind"] == "const":
            e[avail] = mdl["value"]
        else:
            e[avail] = np.interp(sig[avail].to_numpy(float),
                                 mdl["x"], mdl["y"],
                                 left=mdl["y"][0], right=mdl["y"][-1])
        # hard overrides: failures a sigma-only model cannot see
        if k == "EQ":
            e[db[f"sig_{k}"] > EQ_SENTINEL] = np.nan          # source disclaims it
        if k == "DGPMS":
            e[db["is_quinone"] == 1] = np.nan                  # section 2
        out[f"ehat_{k}"] = e
    return out


# ---------------------------------------------------------------- assignment
def assign(db: pd.DataFrame, ehat: pd.DataFrame, tol: float) -> pd.DataFrame:
    cols = [f"ehat_{k}" for k in SOURCES.values()]
    E = ehat[cols].to_numpy(float)
    keys = [c.replace("ehat_", "") for c in cols]
    with np.errstate(invalid="ignore"):
        best = np.nanargmin(np.where(np.isnan(E), np.inf, E), axis=1)
    bestval = np.where(np.all(np.isnan(E), axis=1), np.nan,
                       np.nanmin(np.where(np.isnan(E), np.inf, E), axis=1))
    chosen = np.array([keys[b] for b in best], dtype=object)
    chosen[~np.isfinite(bestval)] = None
    keep = np.isfinite(bestval) & (bestval <= tol)
    dg = np.array([db.iloc[i].get(f"dg_{chosen[i]}", np.nan) if chosen[i] else np.nan
                   for i in range(len(db))], dtype=float)
    op = np.array([db.iloc[i].get(f"op_{chosen[i]}", None) if chosen[i] else None
                   for i in range(len(db))], dtype=object)
    return pd.DataFrame({"chosen_source": np.where(keep, chosen, None),
                         "merged_dg": np.where(keep, dg, np.nan),
                         "merged_operator": np.where(keep, op, None),
                         "ehat": bestval, "kept": keep}, index=db.index)


# ---------------------------------------------------------------- baselines
def baseline_priority(db: pd.DataFrame) -> np.ndarray:
    """dev's Promote_Reaction_Thermodynamics_to_Canonical.py policy: the
    mechanistic tier first, lowest reported error within a tier."""
    out = np.full(len(db), None, dtype=object)
    for tier in (["EQ", "GC"], ["DGPMS"]):
        avail = {k: db[f"dg_{k}"].notna().to_numpy() for k in tier}
        sig = {k: db[f"sig_{k}"].to_numpy(float) for k in tier}
        for i in range(len(db)):
            if out[i] is not None:
                continue
            cand = [k for k in tier if avail[k][i]
                    and not (k == "EQ" and sig[k][i] > EQ_SENTINEL)]
            if cand:
                out[i] = min(cand, key=lambda k: sig[k][i])
    return out


def realized_error(db: pd.DataFrame, chosen: np.ndarray) -> np.ndarray:
    e = np.full(len(db), np.nan)
    for i in range(len(db)):
        k = chosen[i]
        if k and np.isfinite(db.iloc[i].get(f"dg_{k}", np.nan)):
            e[i] = abs(db.iloc[i][f"dg_{k}"] - db.iloc[i]["tecrdb"])
    return e


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    db = load_db()
    print(f"non-EMPTY reactions: {len(db)}")
    for k in SOURCES.values():
        print(f"  {k:6s} covers {db[f'dg_{k}'].notna().sum():6d}")
    union = db[[f"dg_{k}" for k in SOURCES.values()]].notna().any(axis=1)
    print(f"  union: {union.sum()} ({union.mean():.1%}) -- the coverage ceiling")

    truth = load_truth(db)
    print(f"\nTECRDB stereo-exact ground truth joined: {len(truth)} reactions")

    # ---- validation: fit on train, score held-out against baselines
    idx = RNG.permutation(len(truth))
    cut = int(len(truth) * 0.7)
    tr, te = truth.iloc[idx[:cut]].copy(), truth.iloc[idx[cut:]].copy()
    db_tr = db[~db.rxn.isin(te.rxn)]          # keep held-out rxns out of silver too
    models = fit_error_models(tr, db_tr)
    print("\ncalibration (gold = TECRDB train split; silver = vs a trusted-sigma source):")
    for k, m in models.items():
        if m["kind"] == "isotonic":
            print(f"  {k:6s} gold n={m['n_gold']:4d} (median|err| {m['gold_median_err']:.2f}, "
                  f"sigma p90 {m['gold_sigma_p90']:.2f})  silver n={m['n_silver']:6d}  "
                  f"rho(sigma,|err|)={m['spearman_sigma_vs_err']:+.3f}")
        else:
            print(f"  {k:6s} constant fallback {m['value']:.2f}")

    te_e = predict_error(te, models)
    asg = assign(te, te_e, tol=np.inf)          # score the CHOICE, not the filter
    strategies = {
        "assignment (this script)": asg["chosen_source"].to_numpy(),
        "always eQuilibrator": np.where(te["dg_EQ"].notna()
                                        & (te["sig_EQ"] <= EQ_SENTINEL), "EQ", None),
        "always dGPredictor-MS": np.where(te["dg_DGPMS"].notna(), "DGPMS", None),
        "always Group contribution": np.where(te["dg_GC"].notna(), "GC", None),
        "dev priority (EQ>GC, then ML)": baseline_priority(te),
    }
    print(f"\nheld-out TECRDB validation (n={len(te)}), |chosen source - experiment|:")
    rows = []
    for lab, ch in strategies.items():
        err = realized_error(te, ch)
        ok = np.isfinite(err)
        rows.append({"strategy": lab, "n_scored": int(ok.sum()),
                     "median_abs_err": float(np.median(err[ok])),
                     "mean_abs_err": float(np.mean(err[ok])),
                     "p90_abs_err": float(np.percentile(err[ok], 90))})
        print(f"  {lab:32s} n={ok.sum():4d}  median={np.median(err[ok]):5.2f}  "
              f"mean={np.mean(err[ok]):6.2f}  p90={np.percentile(err[ok], 90):6.2f}")
    val = pd.DataFrame(rows)

    # ---- refit on all ground truth, apply database-wide
    models = fit_error_models(truth, db)
    ehat = predict_error(db, models)
    print("\n=== frontier: coverage vs expected-error tolerance ===")
    fr = []
    for tol in (0.5, 1, 1.5, 2, 3, 5, 7.5, 10, 20, 1e9):
        a = assign(db, ehat, tol)
        mix = a.loc[a.kept, "chosen_source"].value_counts()
        fr.append({"tolerance": tol, "n": int(a.kept.sum()),
                   "coverage_of_db": a.kept.mean(),
                   **{f"n_{k}": int(mix.get(k, 0)) for k in SOURCES.values()}})
        lab = "no limit" if tol > 1e8 else f"<= {tol:g}"
        print(f"  ehat {lab:9s} n={a.kept.sum():6d} ({a.kept.mean():5.1%} of db)   "
              + "  ".join(f"{k}:{mix.get(k, 0):6d}" for k in SOURCES.values()))
    pd.DataFrame(fr).to_csv(OUT / "source_assignment_frontier.tsv", sep="\t",
                            index=False, float_format="%.5f")

    a = assign(db, ehat, TOLERANCE)
    out = pd.concat([db[["rxn", "name", "ec", "status", "is_quinone"]],
                     ehat, a], axis=1)
    out.to_csv(ASSIGN_TSV, sep="\t", index=False, float_format="%.4f")
    MODELS_JSON.write_text(json.dumps({
        "tolerance": TOLERANCE, "eq_sentinel": EQ_SENTINEL,
        "models": models, "validation": val.to_dict("records"),
        "note": "ehat is expected |error vs TECRDB| from each source's own sigma; "
                "hard overrides: eQ sentinels and dGPredictor on quinones are never assigned.",
    }, indent=1))
    print(f"\n=== shipped: expected error <= {TOLERANCE} kcal/mol ===")
    mix = a.loc[a.kept, "chosen_source"].value_counts()
    print(f"  {a.kept.sum()} reactions ({a.kept.mean():.1%} of the database, "
          f"{a.kept.sum()/union.sum():.1%} of the union ceiling)")
    for k, v in mix.items():
        print(f"    {k:6s} {v:6d}")
    print(f"\nwrote {ASSIGN_TSV.name}, source_assignment_frontier.tsv, {MODELS_JSON.name}")


def load_assignment(path: Path | None = None) -> pd.DataFrame:
    """rxn -> chosen source, merged dG and operator. Consumer entry point."""
    path = path or ASSIGN_TSV
    if not path.exists():
        print(f"  WARNING: {path} missing; run optimize_thermo_source_assignment.py")
        return pd.DataFrame(columns=["rxn", "chosen_source", "merged_dg",
                                     "merged_operator", "ehat", "kept"])
    return pd.read_csv(path, sep="\t", low_memory=False)


if __name__ == "__main__":
    main()
