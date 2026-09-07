"""Generate Gibbs energies for every ModelSEED reaction and compound.

Production run against the cache and parameters selected in CACHE_STRATEGY.md:
`data/modelseed_cache_b/compounds.sqlite` + `data/cc_params_dedup.npz`.

Every reaction in the database gets a row, including the ones we cannot
compute -- a reaction absent from the output is indistinguishable from one that
failed, and the failure reasons are the map of what curation would buy.

Conditions match every existing ModelSEED thermodynamics script: pH 7,
ionic strength 0.25 M, pMg 3.0, 298.15 K. They are recorded in the output
header so a file can never be read under the wrong convention.
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

import re as _re
# ModelSEED stores all thermodynamics in kcal/mol
KCAL = 4.184
_MISSING = _re.compile(r"(seed:cpd\d+) was not found in the compound cache")


def load_reactions():
    """Every ModelSEED reaction, with a reason when it cannot be posed."""
    out = []
    for f in sorted(glob.glob(f"{MSD}/reaction_*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            rid = r.get("id")
            if not rid:
                continue
            name = (r.get("name") or "").strip()
            stoich = r.get("stoichiometry") or ""
            if not stoich.strip():
                out.append((rid, name, None, "no stoichiometry"))
                continue
            # Sum each compound's coefficient ACROSS compartments. For a
            # transport reaction this leaves only the net chemistry: a pure
            # translocation cancels to nothing (chemical dG = 0), while a
            # driven pump keeps the ATP hydrolysis that powers it. The
            # electrochemical work term depends on the membrane potential and
            # concentrations, so it is not part of a standard dG and is not
            # computed here.
            net, reason = {}, None
            for term in stoich.split(";"):
                parts = term.split(":")
                if len(parts) < 3:
                    reason = "unparseable stoichiometry"
                    break
                try:
                    coeff = float(parts[0])
                except ValueError:
                    reason = "unparseable stoichiometry"
                    break
                net[parts[1]] = net.get(parts[1], 0.0) + coeff
            left, right = [], []
            for cpd, coeff in net.items():
                if abs(coeff) < 1e-9:
                    continue
                side = right if coeff > 0 else left
                n = abs(coeff)
                side.append(f"{n:g} seed:{cpd}" if n != 1 else f"seed:{cpd}")
            if not reason and not left and not right:
                reason = "translocation only (no net chemistry)"
            if reason:
                out.append((rid, name, None, reason))
            elif not (left and right):
                reason = "one-sided reaction"
                out.append((rid, name, None, reason))
            else:
                out.append((rid, name, " + ".join(left) + " = " + " + ".join(right), None))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path,
                    default=ROOT / "data/modelseed_cache_b/compounds.sqlite")
    ap.add_argument("--params", default=str(ROOT / "data/cc_params_dedup.npz"),
                    help='path to a cc_params npz, or "none" to use the '
                         'parameters eQuilibrator ships with')
    ap.add_argument("--out", type=Path, default=ROOT / "data/modelseed_energies.tsv")
    ap.add_argument("--compounds-out", type=Path,
                    default=ROOT / "data/modelseed_formation_energies.tsv")
    ap.add_argument("--p-h", type=float, default=7.0)
    ap.add_argument("--ionic-strength", type=float, default=0.25)
    # eQuilibrator's own default, and what every ModelSEED thermodynamics
    # script uses (they set p_h/ionic_strength/temperature and leave p_mg
    # alone). pMg 3.0 is ~1 mM free Mg2+, roughly physiological. Do not
    # confuse this with the p_mg=14 the *training data* loader substitutes for
    # measurements that never recorded magnesium -- that is missing metadata,
    # not a prediction condition.
    ap.add_argument("--p-mg", type=float, default=3.0)
    ap.add_argument("--temperature", type=float, default=298.15)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    from component_contribution.parameters import CCModelParameters
    from component_contribution.predict import GibbsEnergyPredictor
    from equilibrator_api import ComponentContribution, Q_
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
    print(f"RMSE_inf = {RMSE_INF:,.0f} kJ/mol "
          f"({RMSE_INF / KCAL:,.0f} kcal/mol) -- the scale eQuilibrator uses "
          f"to mean 'unconstrained'", flush=True)

    reactions = load_reactions()
    if args.limit:
        reactions = reactions[:args.limit]
    print(f"ModelSEED reactions: {len(reactions)}", flush=True)

    header = (f"# cache={args.cache} params={args.params} "
              f"p_h={args.p_h} ionic_strength={args.ionic_strength}M "
              f"p_mg={args.p_mg} T={args.temperature}K")
    stats = Counter()
    missing = Counter()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        fh.write(header + "\n")
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reaction_id", "name", "status", "dg_prime_kJ_per_mol",
                    "uncertainty_kJ_per_mol", "dg_prime_kcal_per_mol",
                    "uncertainty_kcal_per_mol", "formula"])
        for i, (rid, name, formula, skip) in enumerate(reactions, 1):
            if i % 5000 == 0:
                print(f"  {i}/{len(reactions)}", flush=True)
            if skip:
                stats[skip] += 1
                w.writerow([rid, name, skip, "", "", "", "", ""])
                continue
            try:
                rxn = cc.parse_reaction_formula(formula)
            except Exception as exc:
                # parse_reaction_formula RAISES for an unknown accession rather
                # than returning a None entry, so a `None in rxn.sparse` check
                # never fires -- every missing compound would otherwise be
                # filed as a parse error. Recover the accession from the
                # message: these are the compounds curation would need to add.
                m = _MISSING.search(str(exc))
                if m:
                    acc = m.group(1)
                    missing[acc] += 1
                    stats["compound not in cache"] += 1
                    w.writerow([rid, name, f"compound not in cache: {acc}",
                                "", "", "", "", formula])
                else:
                    stats["unparseable"] += 1
                    w.writerow([rid, name, "unparseable", "", "", "", "", formula])
                continue
            if not rxn.is_balanced(ignore_atoms=("H", "e-")):
                stats["unbalanced"] += 1
                w.writerow([rid, name, "unbalanced", "", "", "", "", formula])
                continue
            # Ask the predictor WHICH branch it would take before taking it.
            #
            # GibbsEnergyPredictor.standard_dg does:
            #
            #     mu, sigma_fin, sigma_inf, residual = get_reaction_prediction(rxn)
            #     if residual:
            #         return Q_(0, "kJ/mol").plus_minus(RMSE_inf)   # mu DISCARDED
            #     else:
            #         std = |sigma_fin| + RMSE_inf * |sigma_inf|
            #         return Q_(mu, "kJ/mol").plus_minus(max(RMSE_eps, std))
            #
            # Those are two different outcomes and only one of them is a
            # failure. When `residual` is non-empty some reactant cannot be
            # expressed in the model's basis at all, mu is thrown away, and the
            # returned value is literally zero -- plus, after
            # standard_dg_prime, the pH/ionic-strength transform. Storing that
            # is storing the transform and calling it an energy, which is the
            # defect this whole regeneration exists to remove.
            #
            # When `residual` is empty the prediction is real. The uncertainty
            # can still be enormous, because RMSE_inf (1e5 kJ/mol) multiplies
            # any component along an unconstrained direction -- but a real
            # estimate with an honest huge error bar is a result, not a
            # failure, and withholding it is our judgement substituting for the
            # user's.
            #
            # An earlier version keyed on `err > 1e4`, which cannot tell the
            # two apart: both come back with a vast sigma. It therefore
            # discarded ~2,600 genuine estimates along with the empty ones.
            #
            # NOTE the input. ComponentContribution.standard_dg_prime splits
            # the reaction FIRST:
            #
            #     residual_reaction, stored_dg_prime = \
            #         reaction.separate_stored_dg_prime(...)
            #     return predictor.standard_dg_prime(residual_reaction, ...) \
            #            + stored_dg_prime
            #
            # so the predictor never sees compounds that carry a measured dG.
            # Asking get_reaction_prediction about the FULL reaction reports
            # those measured compounds as undecomposable and declines a
            # reaction the library would happily have answered. Split the same
            # way before asking.
            try:
                _residual_rxn, _stored = rxn.separate_stored_dg_prime(
                    p_h=cc.p_h, p_mg=cc.p_mg,
                    ionic_strength=cc.ionic_strength,
                    temperature=cc.temperature)
                _mu, _s_fin, _s_inf, residual = \
                    cc.predictor.get_reaction_prediction(_residual_rxn)
            except Exception:
                stats["prediction failed"] += 1
                w.writerow([rid, name, "prediction failed", "", "", "", "", formula])
                continue

            if residual:
                # No estimate exists. Name the culprits so curation can act.
                culprits = ",".join(sorted(
                    str(getattr(c, "id", c)) for c in residual)) or "unknown"
                stats["undecomposable"] += 1
                w.writerow([rid, name, f"undecomposable: {culprits}",
                            "", "", "", "", formula])
                continue

            try:
                dg = cc.standard_dg_prime(rxn)
            except Exception:
                stats["prediction failed"] += 1
                w.writerow([rid, name, "prediction failed", "", "", "", "", formula])
                continue

            mu, err = dg.value.magnitude, dg.error.magnitude
            stats["ok"] += 1
            # Record the TRUE sigma rather than the string "inf". RMSE_inf is
            # 1e5 kJ/mol, so a reaction with an unconstrained direction comes
            # back at ~that scale; writing "inf" threw away the difference
            # between "unconstrained" and "merely imprecise" and made the two
            # impossible to tell apart downstream. Counted here so the split is
            # visible in the run summary.
            if err >= 0.5 * RMSE_INF:
                stats["  ...of which effectively unconstrained"] += 1
            elif err > 100.0:
                stats["  ...of which very low confidence"] += 1
            w.writerow([rid, name, "ok", f"{mu:.3f}", f"{err:.3f}",
                        f"{mu / KCAL:.3f}", f"{err / KCAL:.3f}", formula])

    if missing:
        # derive from the output name, or a second run with a different --out
        # silently clobbers the first run's list
        mp = args.out.with_name(args.out.stem + "_missing_compounds.tsv")
        with mp.open("w", newline="") as fh:
            mw = csv.writer(fh, delimiter="\t")
            mw.writerow(["seed_id", "reactions_blocked"])
            for acc, n in missing.most_common():
                mw.writerow([acc.replace("seed:", ""), n])
        print(f"wrote {mp}  ({len(missing)} compounds blocking "
              f"{sum(missing.values())} reactions)")

    print(f"\nwrote {args.out}")
    tot = sum(stats.values())
    for k, v in stats.most_common():
        print(f"  {k:<28} {v:6d}  ({100*v/tot:4.1f}%)")


if __name__ == "__main__":
    main()
