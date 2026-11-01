#!/usr/bin/env python
"""Ferredoxin at one electron per unit: cpd11620 (reduced) charge 4 -> 5, the reactions rebalanced.

cpd11620/cpd11621 (Reducedferredoxin/Oxidizedferredoxin, both Fe2R4S6) were charged 4/6, a
2-electron unit, while the database's second ferredoxin pair cpd28082/cpd27757 and the ModelSEED
model templates' compcompounds carry 5/6, one electron per [2Fe-2S]/[4Fe-4S] cluster -- as each
cluster carries -- with Rhea, KEGG, MetaCyc and BiGG writing 2 Fd per 2-electron step. Reactions
written for one convention sit beside reactions written for the other, and in gapfilled models a
pair of them creates electrons (rxn05939 with 1 Fd against rxn14048 with 2: 22 of 403 members of a
methanogenic community exported H2 made from H+ alone).

For every reaction containing cpd11620 or cpd11621 (obsolete ones included), with this
repository's checker (BiochemPy Reactions.balanceReaction):
  - status under 4/6 (asserted equal to the stored status, CK aside), and under 5/6;
  - a reaction OK under 4/6 but not under 5/6 was written in the 2-electron convention: its Fd
    coefficients (cpd11620 and cpd11621 alike) are doubled, which keeps elements balanced (both
    are Fe2R4S6) and restores the charge balance -- applied only if balanceReaction then says OK;
  - every other reaction keeps its stoichiometry and takes its recomputed status (CK kept).
Every changed reaction is re-verified with an independent element and charge counter. Nothing
else changes.

Writes (relative to Biochemistry/):
  compound_*.json/.tsv, reaction_*.json/.tsv           via Compounds/Reactions save methods
  Curation/history/freiburgermsu/cpd_ferredoxin_charge.tsv        the compound edit
  Curation/history/freiburgermsu/rxn_ferredoxin_one_electron.tsv  the coefficient edits
  Curation/history/freiburgermsu/ferredoxin_one_electron_report.tsv  every reaction, before/after
  Curation/history/all_updates.log                                 one line per edit
Run from this directory: ./Ferredoxin_One_Electron.py [--dry-run]
"""
import argparse
import os
import re
import sys
import time
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(HERE, "../../../Libs/Python"))
from BiochemPy import Reactions, Compounds  # noqa: E402

FD_RED, FD_OX = "cpd11620", "cpd11621"
OLD_CHARGE, NEW_CHARGE = 4, 5
CURATOR = "freiburgermsu"


def base(status):
    return "|".join(p for p in status.split("|") if p != "CK")


def independent_balance(stoich, charges, formulas):
    """Element and charge sums by a counter that shares nothing with balanceReaction."""
    el, q = Counter(), 0.0
    for r in stoich:
        c = float(r["coefficient"])
        q += c * charges[r["compound"]]
        for sym, n in re.findall(r"([A-Z][a-z]?)(\d*)", formulas[r["compound"]] or ""):
            if sym != "H":
                el[sym] += c * (int(n) if n else 1)
    return {k: v for k, v in el.items() if abs(v) > 1e-9}, round(q, 6)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    root = os.path.normpath(os.path.join(HERE, "../../../Biochemistry"))
    hist = os.path.join(root, "Curation", "history", CURATOR)

    ch = Compounds()
    cpds = ch.loadCompounds()
    assert cpds[FD_RED]["charge"] == OLD_CHARGE, cpds[FD_RED]["charge"]
    assert cpds[FD_RED]["formula"] == cpds[FD_OX]["formula"], (cpds[FD_RED]["formula"], cpds[FD_OX]["formula"])
    assert cpds[FD_OX]["charge"] == 6

    rh = Reactions()
    rxns = rh.loadReactions()
    charges = {k: v["charge"] for k, v in cpds.items()}
    formulas = {k: v["formula"] for k, v in cpds.items()}
    charges_new = dict(charges, **{FD_RED: NEW_CHARGE})

    rows, coef_edits, outcome = [], [], Counter()
    for rid in sorted(rxns):
        r = rxns[rid]
        st = r["stoichiometry"]
        if not any(x["compound"] in (FD_RED, FD_OX) for x in st):
            continue
        stored = r["status"]
        old = rh.balanceReaction([dict(x) for x in st])
        assert base(old) == base(stored), f"{rid}: checker says {old}, stored {stored}"
        st5 = [dict(x, charge=NEW_CHARGE) if x["compound"] == FD_RED else dict(x) for x in st]
        new = rh.balanceReaction([dict(x) for x in st5])
        doubled = False
        if base(old) == "OK" and new != "OK":
            st_d = [dict(x, coefficient=2 * x["coefficient"]) if x["compound"] in (FD_RED, FD_OX) else dict(x)
                    for x in st5]
            if rh.balanceReaction([dict(x) for x in st_d]) == "OK":
                el, q = independent_balance(st_d, charges_new, formulas)
                assert not el and q == 0, f"{rid}: independent check {el} {q}"
                st5, new, doubled = st_d, "OK", True
                for x in st_d:
                    if x["compound"] in (FD_RED, FD_OX):
                        coef_edits.append((rid, x["compound"], x["coefficient"], x["compartment"]))
        kind = ("doubled Fd" if doubled else f"{base(old)} -> {new}" if base(old) != new else f"unchanged {new}")
        outcome[("obsolete " if r.get("is_obsolete") else "") + (
            "doubled Fd (OK -> OK)" if doubled else
            "newly OK" if new == "OK" and base(old) != "OK" else
            "OK -> not OK" if base(old) == "OK" and new != "OK" else
            "still OK" if new == "OK" else "imbalanced either way")] += 1
        rows.append((rid, int(bool(r.get("is_obsolete"))), stored, Reactions.preserveCK(new, stored), int(doubled),
                     r["equation"]))
        r["stoichiometry"] = st5
        r["status"] = Reactions.preserveCK(new, stored)
        rh.rebuildReaction(r)

    for k, v in sorted(outcome.items()):
        print(f"{v:5d}  {k}")
    print(f"{len(rows)} reactions contain the pair; {len(coef_edits) // 2} Fd coefficient pairs doubled")
    if a.dry_run:
        return
    cpds[FD_RED]["charge"] = NEW_CHARGE
    ch.saveCompounds(cpds)
    rh.saveReactions(rxns)
    os.makedirs(hist, exist_ok=True)
    with open(os.path.join(hist, "cpd_ferredoxin_charge.tsv"), "w") as fh:
        fh.write("ID\tCHARGE\n" + f"{FD_RED}\t{NEW_CHARGE}\n")
    with open(os.path.join(hist, "rxn_ferredoxin_one_electron.tsv"), "w") as fh:
        fh.write("reaction\tcompound\tstoichiometry\tcompartment\n")
        for e in coef_edits:
            fh.write("\t".join(map(str, e)) + "\n")
    with open(os.path.join(hist, "ferredoxin_one_electron_report.tsv"), "w") as fh:
        fh.write("id\tis_obsolete\tstatus_before\tstatus_after\tfd_doubled\tequation_before\n")
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")
    day = time.strftime("%Y-%m-%d", time.gmtime())
    with open(os.path.join(root, "Curation", "history", "all_updates.log"), "a") as fh:
        fh.write("\t".join([day, CURATOR, FD_RED, "charge", str(OLD_CHARGE), str(NEW_CHARGE)]) + "\n")
        for rid, cpd, coef, comp in coef_edits:
            fh.write("\t".join([day, CURATOR, rid, f"{cpd}[{comp}] coefficient", str(coef // 2), str(coef)]) + "\n")


if __name__ == "__main__":
    main()
