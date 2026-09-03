#!/usr/bin/env python
"""Gold / silver / bronze evidence grading on the PROTONATION axis.

This grades one dimension only -- how well founded each reaction's protonation
model is. A full evidence grade would also weigh the energy source (measured
TECRDB anchor vs fitted prediction), the reported uncertainty, and whether the
sources agree on direction. Those are separate axes and are not combined here;
presenting this as the whole grade would overstate it.

Compound grades, in descending order of evidence:

  GOLD    a measured macroscopic ladder (Alberty, or IUPAC digitized from
          Serjeant/Perrin), OR the compound has no ionizable groups, in which
          case there is no protonation model to get wrong.
  SILVER  a ChemAxon macroscopic ladder carried in the pinned eQuilibrator
          cache. A prediction, but a sequential one: it reproduces measured
          literature ladders about 79% of the time, against MolGpKa's 35%.
  BRONZE  MolGpKa microscopic per-site values consumed as a ladder. This is a
          category error independent of accuracy -- eQuilibrator slices the
          list as a sequential ladder, which per-site values are not.
  BRONZE-COLLAPSED
          as bronze, and the compound additionally shows the symmetry-collapse
          signature, so its second and subsequent dissociations are known to be
          copies of the first rather than estimates of anything.

Reaction grade is the MINIMUM over its compounds, which also captures the
mixed-enumeration defect: a reaction with one cache-enumerated and one
MolGpKa-enumerated partner grades bronze, and that inconsistency is exactly what
moved rxn08789 by 355 kcal/mol.

Grades describe the CURRENT build (`cache_final`). The three-tier structural
match is not yet implemented, so a projected column reports what each reaction
would grade once it is.
"""
import csv, pickle, re, sqlite3, sys
from collections import Counter, defaultdict
from pathlib import Path

from paths import (EQ, ZENODO_CACHE as ZEN, THERMO_CACHE, DATA_OUT as OUT,
                   RANKED, REACTIONS, RESOLVED, PRIORITY, require)

CPD = re.compile(r"cpd\d{5}")
ORDER = ["GOLD", "SILVER", "BRONZE", "BRONZE-COLLAPSED"]
RANK = {g: i for i, g in enumerate(ORDER)}

sys.path.insert(0, str(EQ / "tools"))
from modelseed_pkas import cache_seed_identifiers  # noqa: E402


def main():
    fin = str(require(THERMO_CACHE, "rebuilt ModelSEED cache"))
    sm = cache_seed_identifiers(fin)
    cf = sqlite3.connect(fin)
    cz = sqlite3.connect(str(require(ZEN, "upstream Zenodo cache")))

    zl = {}
    for ik, dc in cz.execute(
            "select inchi_key,dissociation_constants from compounds "
            "where inchi_key is not null"):
        v = pickle.loads(dc) if isinstance(dc, (bytes, bytearray)) else dc
        zl[ik] = tuple(round(x, 4) for x in sorted(v or [], reverse=True))

    tool = {}
    with require(RESOLVED, "resolved pKa table").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tool[r["external_id"]] = r["tool"]

    collapsed = set()
    with require(RANKED, "polyprotic ranking").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            collapsed.add(r["seed_id"])

    # Structural-match classification from check_mobile_h.py: which compounds
    # have a Zenodo ladder reachable by a safe match, and which are refused on
    # a fixed-H (tautomer) mismatch.
    rejected, recoverable = set(), set()
    cls = OUT / "structural_match_classification.tsv"
    if cls.exists():
        with cls.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                (recoverable if r["safe"] == "1" else rejected).add(r["seed_id"])
    else:
        print("WARNING: no structural_match_classification.tsv; "
              "run check_mobile_h.py first. Projection will be empty.",
              file=sys.stderr)

    grade, why = {}, {}
    for sid, (cid, ik) in sm.items():
        row = cf.execute(
            "select dissociation_constants from compounds where id=?", (cid,)).fetchone()
        if not row:
            continue
        v = pickle.loads(row[0]) if isinstance(row[0], (bytes, bytearray)) else row[0]
        stored = tuple(round(x, 4) for x in sorted(v or [], reverse=True))
        t = tool.get(sid)
        if not stored:
            grade[sid], why[sid] = "GOLD", "no ionizable groups"
        elif zl.get(ik) == stored:
            grade[sid], why[sid] = "SILVER", "ChemAxon ladder from pinned cache"
        elif t in ("alberty", "iupac"):
            grade[sid], why[sid] = "GOLD", f"measured macroscopic ladder ({t})"
        elif sid in rejected:
            grade[sid], why[sid] = "BRONZE-COLLAPSED", "tautomer match refused; MolGpKa fallback"
        elif sid in collapsed:
            grade[sid], why[sid] = "BRONZE-COLLAPSED", "MolGpKa with symmetry collapse"
        else:
            grade[sid], why[sid] = "BRONZE", "MolGpKa, sites chemically distinct"

    # Projected: under the three-tier structural match, every compound in the
    # recoverable set gains the cache's ChemAxon ladder, i.e. silver. Refused
    # compounds keep the MolGpKa fallback and stay where they are.
    projected = dict(grade)
    for sid in recoverable:
        if sid in projected and RANK[projected[sid]] > RANK["SILVER"]:
            projected[sid] = "SILVER"

    # ---- second scheme: graded by whether the uncertainty can MOVE the answer.
    # Provenance alone puts 79.5% of reactions in bronze, which does not
    # discriminate. FINDINGS measured that a pKa differing near neutrality moves
    # dG'0 about 8x more than one far from it, because a pKa far from the
    # reported pH leaves the microspecies distribution saturated. So grade on
    # whether a compound's UNCERTAIN pKas can sit in the pH 6-8 window.
    #
    # A collapsed compound counts as near-neutral WHATEVER its predicted values,
    # because collapse means the real dissociations are unknown -- MolGpKa
    # returns 2.11 for orthophosphate, whose true pKa2 is 7.20. Skipping that
    # correction grades 2,843 compounds and 10,233 reactions too generously.
    mech = {}
    for sid, (cid, ik) in sm.items():
        row = cf.execute(
            "select dissociation_constants from compounds where id=?", (cid,)).fetchone()
        if not row:
            continue
        v = pickle.loads(row[0]) if isinstance(row[0], (bytes, bytearray)) else row[0]
        st = tuple(round(x, 4) for x in sorted(v or [], reverse=True))
        macro = grade.get(sid) in ("GOLD", "SILVER") and sid not in rejected
        if macro:
            mech[sid] = "MACRO"
        elif any(6.0 <= x <= 8.0 for x in st) or sid in collapsed:
            mech[sid] = "UNCERTAIN_NEAR"
        else:
            mech[sid] = "UNCERTAIN_FAR"

    priority = set()
    pf = PRIORITY
    if pf.exists():
        priority = {l.strip() for l in pf.read_text().split("\n") if l.strip()}

    rows = []
    with require(REACTIONS, "regenerated reaction energies").open() as fh:
        rd = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        for r in rd:
            if r.get("status") != "ok":
                continue
            cps = sorted(set(CPD.findall(r.get("formula") or "")))
            graded = [(c, grade[c]) for c in cps if c in grade]
            if not graded:
                continue
            worst = max(graded, key=lambda x: RANK[x[1]])
            pgraded = [projected.get(c, grade[c]) for c in cps if c in grade]
            pworst = max(pgraded, key=lambda g: RANK[g])
            mc = [mech[c] for c in cps if c in mech]
            mgrade = ("BRONZE" if any(x == "UNCERTAIN_NEAR" for x in mc)
                      else "SILVER" if any(x == "UNCERTAIN_FAR" for x in mc)
                      else "GOLD")
            rows.append({
                "reaction_id": r["reaction_id"], "name": r.get("name", "")[:60],
                "grade": worst[1], "projected_grade": pworst,
                "mechanism_grade": mgrade,
                "limiting_compound": worst[0],
                "reason": why[worst[0]],
                "n_compounds": len(cps), "n_graded": len(graded),
                "distinct_grades": len({g for _, g in graded}),
                "in_priority_scope": int(r["reaction_id"] in priority),
                "dg_prime_kcal_per_mol": r.get("dg_prime_kcal_per_mol", ""),
                "uncertainty_kcal_per_mol": r.get("uncertainty_kcal_per_mol", ""),
            })

    with (OUT / "reaction_evidence_grades.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    tot = Counter(r["grade"] for r in rows)
    pri = Counter(r["grade"] for r in rows if r["in_priority_scope"])
    n, npri = len(rows), sum(1 for r in rows if r["in_priority_scope"])
    print("PROTONATION EVIDENCE GRADE, scored reactions in cache_final\n")
    print(f"{'grade':20s} {'all':>8s} {'':>7s}  {'priority scope':>14s}")
    for g in ORDER:
        print(f"  {g:18s} {tot[g]:>8,} {100*tot[g]/n:6.1f}%  "
              f"{pri[g]:>8,} {100*pri[g]/npri if npri else 0:6.1f}%")
    print(f"  {'TOTAL':18s} {n:>8,}          {npri:>8,}")
    pro = Counter(r["projected_grade"] for r in rows)
    print("\nPROJECTED, after the three-tier structural match:")
    for g in ORDER:
        d = pro[g] - tot[g]
        print(f"  {g:18s} {pro[g]:>8,} {100*pro[g]/n:6.1f}%   "
              f"({d:+,} vs current)")
    moved = sum(1 for r in rows if r["grade"] != r["projected_grade"])
    print(f"  reactions improving a grade: {moved:,} ({100*moved/n:.1f}%)")

    mg = Counter(r["mechanism_grade"] for r in rows)
    mgp = Counter(r["mechanism_grade"] for r in rows if r["in_priority_scope"])
    print("\nMECHANISM SCHEME -- can the uncertainty reach the pH 6-8 window?")
    for g in ("GOLD", "SILVER", "BRONZE"):
        print(f"  {g:18s} {mg[g]:>8,} {100*mg[g]/n:6.1f}%  "
              f"{mgp[g]:>8,} {100*mgp[g]/npri if npri else 0:6.1f}%")

    mixed = sum(1 for r in rows if r["distinct_grades"] > 1)
    print(f"\nreactions mixing evidence grades across their compounds: {mixed:,} "
          f"({100*mixed/n:.1f}%)")
    print("\nmost frequent limiting compounds:")
    lim = Counter((r["limiting_compound"], r["reason"]) for r in rows
                  if r["grade"].startswith("BRONZE"))
    for (c, rsn), k in lim.most_common(10):
        print(f"   {c}  {k:>6,} reactions   {rsn}")
    print(f"\nwrote {OUT / 'reaction_evidence_grades.tsv'}")


if __name__ == "__main__":
    main()
