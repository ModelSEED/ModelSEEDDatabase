#!/usr/bin/env python3
"""Recommend WHICH thermodynamic source to use for a reaction.

Separate from ``grade_thermo_sources.py`` on purpose. A grade answers "how much
should I trust this number"; a recommendation answers "which number should I
use". They are different jobs and the same statistic does not serve both.

    magnitude target   s*(i) = argmin_s ehat_s(i)      uncertainty DOES arbitrate
    direction target   s*(i) = first feasible source   uncertainty does NOT
                               in a fixed priority

THE NEGATIVE RESULT THAT SHAPED THIS, stated up front because it is the whole
design. Held-out direction accuracy on the 802 TECRDB stereo-exact anchors,
20 x 70/30 splits, every strategy at 100% coverage unless noted:

    priority EQ > DGPMS > GC          95.9% +/- 1.1     <-- shipped
    eQuilibrator only                 95.9% +/- 1.2     (98.3% coverage)
    priority + risk veto at 0.20      94.2% +/- 1.1
    priority + risk veto at 0.05      93.8% +/- 1.4
    argmin calibrated tau             93.7% +/- 1.1
    argmin ehat (magnitude rule)      93.7% +/- 1.1
    priority + risk veto at 0.02      93.0% +/- 1.3
    dGPredictor-ModelSEED only        91.6% +/- 1.1
    argmin direction risk             90.9% +/- 1.2
    Group contribution only           85.1% +/- 2.3

Every attempt to arbitrate BETWEEN sources using their reported uncertainty
lost to a fixed priority order, and adding a risk veto on top of priority made
things monotonically worse -- the veto only ever pushes a reaction off the
better source onto a worse one.

WHY argmin-risk FAILS. The risk below is P(this source's own call is overturned
by this source's own uncertainty). That is PRECISION, not ACCURACY. A source
whose dG is confidently wrong -- small tau, far from any cascade breakpoint,
wrong region -- scores risk ~ 0 and beats a source that is right but sits near
a band edge. Centering the integral on a biased point estimate cannot detect
the bias. The quantity is well defined and correctly computed; it is simply not
the quantity that selects a source.

WHAT THE UNCERTAINTY IS STILL FOR. Three jobs it does do, all validated:
  * FEASIBILITY -- eQuilibrator sentinels (sigma > 100, the source disclaiming
    the reaction), the MetaNetX collision list, dGPredictor-ModelSEED on
    quinones. These remove a source outright.
  * ABSTENTION -- the risk on the CHOSEN source decides whether to answer at
    all. Within a source it is informative: eQuilibrator's accuracy by its own
    sigma quartile runs 100.0 / 98.6 / 91.9 / 92.0%.
  * MAGNITUDE ARBITRATION -- for the magnitude target argmin ehat is the right
    rule and does beat always-eQuilibrator (mean |error| 1.03 vs 1.75 kcal/mol,
    validated in optimize_thermo_source_assignment.py).

CAVEAT ON THE PRIORITY ORDER. It comes from measured accuracy on TECRDB, but
eQuilibrator's reactant-contribution layer is FITTED on TECRDB and dGPredictor
was trained on 4,001 measurements from it, so both are partly in-sample and
Group Contribution is not. eQuilibrator's edge is concentrated where its own
sigma is smallest (100% in its lowest sigma quartile, tied with
dGPredictor-ModelSEED at 92% in the third), which is what partial memorisation
would look like. Treat "prefer eQuilibrator" as the best rule available on the
evidence we have, not as a settled fact about the methods.

THE DIRECTION RISK (still computed, reported, and used for abstention)
----------------------------------------------------------------------
Risk = 1 - P(the cascade's operator from the TRUE dG equals the operator it
gives from dG_s), computed EXACTLY rather than sampled. The cascade's operator
is a piecewise-constant function of dG; holding stoichiometry fixed it changes
only at these breakpoints (from reversibility_heuristics):

    -dge - RT(pdt_max + rct_min)     stored_bounds max crosses 0
    +dge - RT(pdt_min + rct_max)     stored_bounds min crosses 0
    -2 - RT*rgt_sum, 2 - RT*rgt_sum  the mMdeltaG +/-2 band edges
    -RT*rgt_sum                      mMdeltaG sign flip (low-energy)
    2/points - RT*rgt_sum            low-energy threshold

The REAL cascade is evaluated once per interval, then the calibrated normal
N(dG_s, tau_s^2) is integrated over the intervals whose operator matches -- a
closed-form sum of normal CDFs, exact with respect to the cascade rather than
an approximation of one heuristic. ATP-synthase and ABC-transporter reactions
short-circuit before any dG is read, so their risk is exactly 0.

CALIBRATED SIGMA
----------------
tau_s is NOT the reported sigma -- the three sigma scales are not commensurable
(Group Contribution overstates its error 2.2x, dGPredictor-ModelSEED 1.5x,
eQuilibrator understates 1.6x). tau_s(i) = k_s * ehat_s(i) / sqrt(2/pi): the
calibrated expected |error| converted to a Gaussian scale, times one per-source
scalar k_s fitted so +/-tau actually covers 68.3% of measured errors on the
anchor. Coverage before and after is printed and stored, so the calibration is
checkable rather than asserted.

OUTPUT (results/thermo_recommendation/)
    recommendation_direction.tsv  per reaction: chosen source, dG, risk, why
    recommendation_magnitude.tsv  same, argmin-ehat rule
    recommendation_models.json    tau calibration + the full ablation above
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ANALYSIS_DIR = ANALYSIS_DIR
MSDB_CODE = MSDB_CODE
OUT = Path(os.environ.get("RECOMMEND_OUT",
                          str(ANALYSIS_DIR / "results" / "thermo_recommendation")))
sys.path.insert(0, str(Path(__file__).resolve().parent))
REPO_ROOT = Path(__file__).resolve().parents[3]
MSDB_ROOT = Path(os.environ.get("MSDB_ROOT", str(REPO_ROOT)))
MSDB_CODE = Path(os.environ.get("MSDB_CODE", str(REPO_ROOT)))
ANALYSIS_DIR = Path(os.environ.get("THERMO_GRADING_OUT",
                                   str(REPO_ROOT / "Biochemistry" / "Thermodynamics" / "SourceGrading")))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from reversibility_heuristics import (  # noqa: E402
    DEFAULT_HEURISTICS, RT_CONST, _low_energy_points, _walk_stoichiometry,
    explicit_energy, run_reversibility,
)
from optimize_thermo_source_assignment import (  # noqa: E402
    EQ_SENTINEL, fit_error_models, load_db, predict_error,
)
from grade_thermo_sources import load_tecrdb, load_vetoes  # noqa: E402

K = ["GC", "EQ", "DGPMS"]
LABEL = {"GC": "Group contribution", "EQ": "eQuilibrator",
         "DGPMS": "dGPredictor-ModelSEED"}
SQRT_2_OVER_PI = math.sqrt(2.0 / math.pi)      # E|X| = tau * sqrt(2/pi)
TAU_FLOOR = 0.05                                # kcal/mol
TRUTH_SIGMA = 0.15                              # TECRDB median experimental sd
BIG = 1e6
RNG = np.random.default_rng(20260812)

# Source priority, per target, ordered by MEASURED accuracy on the TECRDB
# anchor -- not by a hunch and not by reported confidence. Direction:
# eQuilibrator 95.5% > dGPredictor-ModelSEED 91.8% > Group Contribution 85.5%.
# Magnitude (median |error|): eQuilibrator 0.45 ~ dGPredictor-ModelSEED 0.47 >
# Group Contribution 1.57 kcal/mol.
PRIORITY_DIRECTION = ["EQ", "DGPMS", "GC"]
PRIORITY_MAGNITUDE = ["EQ", "DGPMS", "GC"]


def _phi(z):
    return 0.5 * (1.0 + np.vectorize(math.erf)(z / math.sqrt(2.0)))


# --------------------------------------------------------------- calibration
def calibrate_tau(db, anchor, ehat) -> dict:
    """Per-source scale k_s so that +/-tau covers the nominal 68.3% of measured
    |error| on the anchor set. Returns {source: {k, coverage_before, after}}."""
    out = {}
    for k in K:
        e = ehat[f"ehat_{k}"]
        raw_tau = np.maximum(e / SQRT_2_OVER_PI, TAU_FLOOR)
        m = anchor.index
        err = (db.loc[m, f"dg_{k}"] - db.loc[m, "tecrdb_dg"]).abs()
        t = raw_tau.loc[m]
        ok = err.notna() & t.notna()
        if ok.sum() < 30:
            out[k] = {"k": 1.0, "n": int(ok.sum()), "coverage_before": None,
                      "coverage_after": None}
            continue
        before = float((err[ok] <= t[ok]).mean())
        # single scalar that lands the empirical coverage on 0.683
        grid = np.geomspace(0.05, 5.0, 400)
        cov = np.array([float((err[ok] <= g * t[ok]).mean()) for g in grid])
        kbest = float(grid[int(np.argmin(np.abs(cov - 0.683)))])
        out[k] = {"k": kbest, "n": int(ok.sum()), "coverage_before": before,
                  "coverage_after": float((err[ok] <= kbest * t[ok]).mean())}
    return out


def tau_for(db, ehat, kcal: dict) -> pd.DataFrame:
    t = pd.DataFrame(index=db.index)
    for k in K:
        t[f"tau_{k}"] = np.maximum(
            kcal[k]["k"] * ehat[f"ehat_{k}"] / SQRT_2_OVER_PI, TAU_FLOOR)
    return t


# ------------------------------------------------- exact direction agreement
def breakpoints(rxn_entry, dge_truth: float):
    """The dG values at which the cascade's operator can change."""
    terms = _walk_stoichiometry(rxn_entry["stoichiometry"])
    a = RT_CONST * (terms["pdt_max"] + terms["rct_min"])
    b = RT_CONST * (terms["pdt_min"] + terms["rct_max"])
    c = RT_CONST * terms["rgt_sum"]
    pts = _low_energy_points(rxn_entry["stoichiometry"], terms["phosphates"])
    bp = [-dge_truth - a, dge_truth - b, -2.0 - c, 2.0 - c, -c]
    if pts:
        bp.append(2.0 / pts - c)
    return sorted(set(round(x, 9) for x in bp if np.isfinite(x)))


def operator_intervals(rxn_entry, dge_truth: float = TRUTH_SIGMA):
    """[(lo, hi, operator), ...] partitioning the dG line.

    The operator comes from the REAL cascade evaluated at an interior point of
    each interval, so this cannot drift from the cascade's behaviour.
    """
    bp = breakpoints(rxn_entry, dge_truth)
    edges = [-BIG] + bp + [BIG]
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi <= lo:
            continue
        mid = 0.5 * (lo + hi)
        _, op, _ = run_reversibility(rxn_entry, explicit_energy(mid, dge_truth),
                                     DEFAULT_HEURISTICS)
        if out and out[-1][2] == op:
            out[-1] = (out[-1][0], hi, op)       # merge equal neighbours
        else:
            out.append((lo, hi, op))
    return out


def p_operator_matches(intervals, op_target: str, mu: float, tau: float) -> float:
    """P(operator(true dG) == op_target) under true dG ~ N(mu, tau^2)."""
    total = 0.0
    for lo, hi, op in intervals:
        if op != op_target:
            continue
        total += float(_phi((hi - mu) / tau) - _phi((lo - mu) / tau))
    return min(max(total, 0.0), 1.0)


# ------------------------------------------------------------------- risks
def direction_risk(reactions, db, tau, ops, feasible) -> pd.DataFrame:
    """1 - P(this source's operator survives its own uncertainty)."""
    risk = pd.DataFrame(np.nan, index=db.index, columns=[f"risk_{k}" for k in K])
    cache_key = None
    for pos, rxn_id in enumerate(db.rxn):
        entry = reactions.get(rxn_id)
        if entry is None or entry.get("status") == "EMPTY":
            continue
        if not any(feasible[k].iat[pos] for k in K):
            continue
        intervals = operator_intervals(entry)
        # ATPS / ABCT short-circuit before any dG is read -> single interval
        certain = len(intervals) == 1
        for k in K:
            if not feasible[k].iat[pos]:
                continue
            op = ops[k].iat[pos]
            if not isinstance(op, str) or not op:
                continue
            if certain:
                risk.iat[pos, K.index(k)] = 0.0
                continue
            mu = float(db[f"dg_{k}"].iat[pos])
            t = float(tau[f"tau_{k}"].iat[pos])
            risk.iat[pos, K.index(k)] = 1.0 - p_operator_matches(intervals, op, mu, t)
        cache_key = rxn_id
    return risk


def source_operators(reactions, db) -> dict:
    """Each source's deterministic cascade operator, per reaction."""
    from reversibility_heuristics import per_source_energy
    ops = {}
    for k in K:
        src = per_source_energy(LABEL[k])
        col = []
        for rxn_id in db.rxn:
            e = reactions.get(rxn_id)
            if e is None or e.get("status") == "EMPTY":
                col.append("")
                continue
            _, op, _ = run_reversibility(e, src, DEFAULT_HEURISTICS)
            col.append(op or "")
        ops[k] = pd.Series(col, index=db.index)
    return ops


def feasibility(db, veto_eq) -> dict:
    f = {}
    for k in K:
        ok = db[f"dg_{k}"].notna() & db[f"sig_{k}"].notna()
        if k == "EQ":
            ok &= db["sig_EQ"] <= EQ_SENTINEL
            ok &= ~db["rxn"].isin(veto_eq)
        if k == "DGPMS":
            ok &= db["is_quinone"] == 0
        f[k] = ok
    return f


# -------------------------------------------------------------------- choose
RISK_EPS = 1e-4        # risks within this are a tie, and very many are


def rank_key(risk_val, tau_val):
    """Sort key for one source on one reaction: lowest risk wins; ties broken
    by the SHARPER calibrated posterior.

    The tie-break is load-bearing, not cosmetic. The direction risk saturates
    at 0 for most reactions -- once a source's dG is far from every cascade
    breakpoint relative to its own tau, the call is certain and every such
    source ties. Risk alone therefore cannot discriminate on the majority of
    the database, and whatever breaks the tie IS the decision rule there.

    Preferring the smaller tau is the honest choice: among sources equally
    unlikely to have their call overturned, take the one whose energy is most
    precisely known on the common calibrated scale. It is decided before
    looking at any outcome, and it is not "prefer eQuilibrator" -- eQ simply
    has the smallest calibrated tau most of the time.
    """
    return (round(float(risk_val) / RISK_EPS), float(tau_val))


def choose_by_priority(db, risk, feasible, tol, priority) -> pd.DataFrame:
    """Pick the highest-priority feasible source; report its risk; abstain on it.

    THIS IS THE RULE THAT SURVIVED VALIDATION FOR THE DIRECTION TARGET, and it
    deliberately does NOT use uncertainty to choose. Every uncertainty-based
    arbitration tested came out worse (see ``validate``): argmin-risk 90.9%,
    argmin-tau 93.7%, argmin-ehat 93.7%, priority 95.9%, all at 100% coverage
    on the TECRDB anchor. Adding a risk veto on top of priority made it worse
    monotonically (94.2% at 0.2, 93.8% at 0.05, 93.0% at 0.02) because the veto
    pushes reactions off the better source onto a worse one.

    Uncertainty still does two jobs here, just not the arbitration job:
      * feasibility -- sentinels, collisions and the quinone veto remove a
        source outright (folded into ``feasible``);
      * abstention  -- ``kept`` is False when even the chosen source's risk
        exceeds the tolerance, so the caller gets nothing rather than a coin
        flip.
    """
    n = len(db)
    chosen = np.full(n, None, dtype=object)
    chosen_risk = np.full(n, np.nan)
    for k in reversed(priority):        # reverse so higher priority overwrites
        m = feasible[k].to_numpy()
        chosen = np.where(m, k, chosen)
        chosen_risk = np.where(m, risk[f"risk_{k}"].to_numpy(), chosen_risk)
    have = chosen != None                                        # noqa: E711
    keep = have & (np.nan_to_num(chosen_risk, nan=1.0) <= tol)
    dg = np.array([db[f"dg_{chosen[i]}"].iat[i] if chosen[i] else np.nan
                   for i in range(n)], dtype=float)
    sig = np.array([db[f"sig_{chosen[i]}"].iat[i] if chosen[i] else np.nan
                    for i in range(n)], dtype=float)
    return pd.DataFrame({
        "chosen_source": np.where(keep, chosen, None),
        "chosen_label": [LABEL[c] if keep[i] and c else None
                         for i, c in enumerate(chosen)],
        "recommended_dg": np.where(keep, dg, np.nan),
        "recommended_sigma": np.where(keep, sig, np.nan),
        "risk": chosen_risk,
        "n_feasible": np.column_stack([feasible[k] for k in K]).sum(axis=1),
        "kept": keep}, index=db.index)


def choose(db, risk, feasible, tol, tau) -> pd.DataFrame:
    R = np.column_stack([np.where(feasible[k], risk[f"risk_{k}"], np.inf) for k in K])
    R = np.where(np.isnan(R), np.inf, R)
    T = np.column_stack([np.where(feasible[k], tau[f"tau_{k}"], np.inf) for k in K])
    T = np.where(np.isnan(T), np.inf, T)
    # lexicographic (risk bucket, tau): quantise risk so near-ties go to tau
    keyed = np.round(R / RISK_EPS) + np.clip(T, 0, 1e9) / 1e12
    keyed = np.where(np.isfinite(R), keyed, np.inf)
    best = np.argmin(keyed, axis=1)
    bestrisk = R[np.arange(len(R)), best]
    chosen = np.array([K[b] for b in best], dtype=object)
    ok = np.isfinite(bestrisk)
    chosen[~ok] = None
    keep = ok & (bestrisk <= tol)
    dg = np.array([db[f"dg_{chosen[i]}"].iat[i] if chosen[i] else np.nan
                   for i in range(len(db))], dtype=float)
    sig = np.array([db[f"sig_{chosen[i]}"].iat[i] if chosen[i] else np.nan
                    for i in range(len(db))], dtype=float)
    return pd.DataFrame({
        "chosen_source": np.where(keep, chosen, None),
        "chosen_label": [LABEL[c] if keep[i] and c else None
                         for i, c in enumerate(chosen)],
        "recommended_dg": np.where(keep, dg, np.nan),
        "recommended_sigma": np.where(keep, sig, np.nan),
        "risk": np.where(ok, bestrisk, np.nan),
        "n_feasible": np.column_stack([feasible[k] for k in K]).sum(axis=1),
        "kept": keep}, index=db.index)


# ---------------------------------------------------------------- validation
def validate(db, anchor_idx, risk, feasible, ops, ref_op, tau, ehat, n_splits=20):
    """Direction accuracy for the recommender and its baselines.

    HELD OUT ONLY IN THE SCORING SENSE. risk, tau and ehat are computed once on
    the full anchor before this function is called, so the 70/30 permutation
    changes which reactions are SCORED, not which are FITTED. The +/- figures
    are therefore test-set sampling variability, not out-of-sample
    generalisation.

    That matters less here than it would elsewhere, for two reasons: the
    strategies being compared are mostly parameter-free (a fixed priority order,
    or argmin over a precomputed quantity), and the only fitted quantities are
    three scalars. But the distinction is real -- see ``validate_cv`` in
    grade_thermo_sources.py, which does refit per fold, and moves the grading
    numbers by at most 0.04 kcal/mol.
    """
    ids = np.array(anchor_idx)
    rows = []
    strategies = ["recommend_direction", "risk_only_no_tau", "tau_only_no_risk",
                  "recommend_magnitude", "eq_only", "dgpms_only", "gc_only",
                  "eq_first_then_risk", "priority_plain",
                  "priority_veto_0.02", "priority_veto_0.05", "priority_veto_0.2"]
    acc = {s: [] for s in strategies}
    cov = {s: [] for s in strategies}
    for _ in range(n_splits):
        perm = RNG.permutation(len(ids))
        te = ids[perm[int(0.7 * len(ids)):]]
        for s in strategies:
            n = ok = 0
            for i in te:
                truth = ref_op.get(db.rxn.iat[i])
                if truth is None:
                    continue
                if s == "recommend_direction":
                    cands = [(rank_key(risk[f"risk_{k}"].iat[i], tau[f"tau_{k}"].iat[i]), k)
                             for k in K if feasible[k].iat[i]
                             and np.isfinite(risk[f"risk_{k}"].iat[i])]
                    if not cands:
                        continue
                    pick = min(cands)[1]
                elif s == "risk_only_no_tau":
                    cands = [(round(float(risk[f"risk_{k}"].iat[i]) / RISK_EPS), k)
                             for k in K if feasible[k].iat[i]
                             and np.isfinite(risk[f"risk_{k}"].iat[i])]
                    if not cands:
                        continue
                    pick = min(cands)[1]
                elif s == "tau_only_no_risk":
                    cands = [(float(tau[f"tau_{k}"].iat[i]), k) for k in K
                             if feasible[k].iat[i] and np.isfinite(tau[f"tau_{k}"].iat[i])]
                    if not cands:
                        continue
                    pick = min(cands)[1]
                elif s == "recommend_magnitude":
                    cands = [(ehat[f"ehat_{k}"].iat[i], k) for k in K
                             if feasible[k].iat[i] and np.isfinite(ehat[f"ehat_{k}"].iat[i])]
                    if not cands:
                        continue
                    pick = min(cands)[1]
                elif s.startswith("priority_veto"):
                    thr = float(s.split("_")[-1])
                    pick = None
                    for cand in PRIORITY_DIRECTION:
                        if not feasible[cand].iat[i]:
                            continue
                        r = risk[f"risk_{cand}"].iat[i]
                        if np.isfinite(r) and r <= thr:
                            pick = cand
                            break
                    if pick is None:
                        rest = [c for c in PRIORITY_DIRECTION if feasible[c].iat[i]]
                        if not rest:
                            continue
                        pick = rest[0]
                elif s == "priority_plain":
                    rest = [c for c in PRIORITY_DIRECTION if feasible[c].iat[i]]
                    if not rest:
                        continue
                    pick = rest[0]
                elif s == "eq_first_then_risk":
                    if feasible["EQ"].iat[i]:
                        pick = "EQ"
                    else:
                        cands = [(rank_key(risk[f"risk_{k}"].iat[i], tau[f"tau_{k}"].iat[i]), k)
                                 for k in K if feasible[k].iat[i]
                                 and np.isfinite(risk[f"risk_{k}"].iat[i])]
                        if not cands:
                            continue
                        pick = min(cands)[1]
                else:
                    pick = {"eq_only": "EQ", "dgpms_only": "DGPMS", "gc_only": "GC"}[s]
                    if not feasible[pick].iat[i]:
                        continue
                op = ops[pick].iat[i]
                if not op:
                    continue
                n += 1
                ok += (op == truth)
            if n:
                acc[s].append(ok / n)
                cov[s].append(n / len(te))
    for s in strategies:
        rows.append({"strategy": s, "mean_accuracy": float(np.mean(acc[s])),
                     "sd_accuracy": float(np.std(acc[s])),
                     "mean_coverage": float(np.mean(cov[s]))})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--target", choices=("direction", "magnitude", "both"),
                    default="both")
    ap.add_argument("--tolerance", type=float, default=None,
                    help="max risk to accept; default 0.10 for direction, "
                         "2.0 kcal/mol for magnitude")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    from build_graded_direction_maps import load_reactions
    db = load_db().reset_index(drop=True)
    tec = db[["rxn"]].merge(load_tecrdb(), on="rxn", how="left")
    tec.index = db.index
    db["tecrdb_dg"] = tec["tecrdb_dg"]
    anchor = db[tec["match_tier"].eq("stereo_exact").to_numpy() & db.tecrdb_dg.notna()]
    print(f"{len(db)} reactions; TECRDB stereo-exact anchor {len(anchor)}")

    reactions = load_reactions()
    veto_eq = load_vetoes()
    feasible = feasibility(db, veto_eq)
    for k in K:
        print(f"  {LABEL[k]:24s} feasible on {int(feasible[k].sum()):6d}")

    ehat = predict_error(db, fit_error_models(
        anchor.rename(columns={"tecrdb_dg": "tecrdb"}), db))
    kcal = calibrate_tau(db, anchor, ehat)
    print("\ntau calibration (target coverage 0.683 of |error| within +/-tau):")
    for k in K:
        c = kcal[k]
        print(f"  {LABEL[k]:24s} k={c['k']:.3f}  coverage {c['coverage_before']:.3f}"
              f" -> {c['coverage_after']:.3f}  (n={c['n']})")
    tau = tau_for(db, ehat, kcal)

    print("\nrunning the cascade per source...")
    ops = source_operators(reactions, db)
    print("computing exact direction risk (interval integration)...")
    risk = direction_risk(reactions, db, tau, ops, feasible)
    for k in K:
        r = risk[f"risk_{k}"].dropna()
        print(f"  {LABEL[k]:24s} n={len(r):6d}  median risk {r.median():.4f}  "
              f"share risk<=0.05 {float((r <= 0.05).mean()):.1%}")

    # reference operator from the experiment, for validation
    from reversibility_heuristics import explicit_energy as _ee
    ref_op = {}
    t2 = load_tecrdb()
    for r in t2[t2.match_tier == "stereo_exact"].itertuples():
        e = reactions.get(r.rxn)
        if e is None or e.get("status") == "EMPTY":
            continue
        _, op, _ = run_reversibility(e, _ee(float(r.tecrdb_dg),
                                            float(r.tecrdb_sd or 0.0)),
                                     DEFAULT_HEURISTICS)
        if op:
            ref_op[r.rxn] = op

    val = validate(db, list(anchor.index), risk, feasible, ops, ref_op, tau, ehat)
    print("\nheld-out direction accuracy on the TECRDB anchor "
          "(20 x 70/30 splits, coverage = share of held-out reactions the "
          "strategy will answer for):")
    for r in val.sort_values("mean_accuracy", ascending=False).itertuples():
        print(f"  {r.strategy:22s} {r.mean_accuracy:6.1%} +/- {r.sd_accuracy:.1%}"
              f"   coverage {r.mean_coverage:5.1%}")

    targets = ["direction", "magnitude"] if args.target == "both" else [args.target]
    summary = {}
    for target in targets:
        if target == "direction":
            rk = risk
            tol = 0.35 if args.tolerance is None else args.tolerance
        else:
            rk = ehat.rename(columns={f"ehat_{k}": f"risk_{k}" for k in K})
            tol = 2.0 if args.tolerance is None else args.tolerance
        # direction: priority selection (validated); magnitude: argmin ehat
        # (validated in optimize_thermo_source_assignment.py -- there the
        # uncertainty DOES arbitrate correctly, mean |error| 1.03 vs 1.75).
        if target == "direction":
            ch = choose_by_priority(db, rk, feasible, tol, PRIORITY_DIRECTION)
        else:
            ch = choose(db, rk, feasible, tol, tau)
        out = pd.concat([db[["rxn", "name", "ec", "status"]],
                         pd.DataFrame({f"dg_{k}": db[f"dg_{k}"] for k in K}),
                         pd.DataFrame({f"risk_{k}": rk[f"risk_{k}"] for k in K}),
                         ch], axis=1)
        # TECRDB overrides: a measurement outranks every prediction
        have = db.tecrdb_dg.notna()
        out.loc[have, "chosen_source"] = "TECRDB"
        out.loc[have, "chosen_label"] = "TECRDB"
        out.loc[have, "recommended_dg"] = db.loc[have, "tecrdb_dg"]
        out.loc[have, "recommended_sigma"] = tec.loc[have, "tecrdb_sd"]
        out.loc[have, "risk"] = 0.0
        out.loc[have, "kept"] = True
        path = OUT / f"recommendation_{target}.tsv"
        out.to_csv(path, sep="\t", index=False, float_format="%.5f")
        mix = out.loc[out.kept, "chosen_label"].value_counts()
        summary[target] = {"tolerance": tol, "n_kept": int(out.kept.sum()),
                           "mix": {k: int(v) for k, v in mix.items()}}
        print(f"\n=== target={target}, tolerance {tol} ===")
        print(f"  {int(out.kept.sum())} reactions recommended")
        for k, v in mix.items():
            print(f"    {k:24s} {v:6d}")
        print(f"  wrote {path.name}")

    (OUT / "recommendation_models.json").write_text(json.dumps({
        "tau_calibration": kcal, "truth_sigma": TRUTH_SIGMA,
        "validation": val.to_dict("records"), "targets": summary}, indent=1))
    print(f"\nwrote {OUT}/recommendation_models.json")


def load_recommendation(target: str = "direction") -> pd.DataFrame:
    path = OUT / f"recommendation_{target}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} missing; run recommend_thermo_source.py")
    return pd.read_csv(path, sep="\t", low_memory=False)


if __name__ == "__main__":
    main()
