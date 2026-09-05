#!/usr/bin/env python
"""Is it safe to match our structures to the eQuilibrator cache modulo /h and /p?

FAD matches Zenodo on formula and connectivity but not on InChIKey, because the
two InChIs partition the mobile hydrogens differently. Ignoring the /h layer
would recover a Marvin ladder for that whole class -- but /h also encodes
tautomers, and a genuine tautomer can have genuinely different pKas. So the
question is whether these are the same protons regrouped, or different protons.

The /h layer splits into a fixed part (atoms carrying non-exchangeable H) and
zero or more mobile-H groups, e.g.::

    ours   h3-4,8-9,...,37-41H,5-7H2,1-2H3(H,44,45)(H,46,47)(H2,28,29,30)(H,34,42,43)
    zenodo h3-4,8-9,...,37-41H,5-7H2,1-2H3(H5,28,29,30,34,42,43,44,45,46,47)

A regrouping is safe iff the fixed part is identical, the total number of mobile
protons is identical, and the SET of atoms they may occupy is identical -- only
the partition into groups differs. Anything else is a real difference in where
the hydrogens are, and the ladder must not be transferred on that basis.

/p is reported but does not disqualify a match: it is a protonation-state
offset, and a ladder spans protonation states by construction. It does change
the reference proton count, which is what the magnesium guard arithmetic uses,
so it is surfaced per compound.
"""

if __name__ == "__main__":
    # Validate arguments BEFORE importing anything or touching the database.
    # These scripts mutate the database, and without this an unknown flag or a
    # mistyped mode was silently ignored and the script ran with its defaults:
    # asking Estimate_Reaction_Reversibility.py for --help rewrote 122 files.
    # Placed above the imports so --help works even where a dependency is
    # missing from the path.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import csv, pickle, re, sqlite3, sys
from collections import defaultdict

from paths import ZENODO_CACHE, STRUCTURES, RANKED, DATA_OUT, require

MOBILE = re.compile(r"\(H(\d*)([-+]?),([\d,]+)\)")


def split_inchi(inchi):
    if not inchi or not inchi.startswith("InChI="):
        return None
    parts = inchi.split("/")
    get = lambda p: next((x for x in parts if x.startswith(p)), "")
    return {"formula": parts[1] if len(parts) > 1 else "",
            "c": get("c"), "h": get("h"), "p": get("p"),
            "t": get("t"), "m": get("m"), "s": get("s")}


def parse_h(h):
    """-> (fixed_spec, total_mobile_protons, frozenset(atoms that may bear them))"""
    if not h:
        return "", 0, frozenset()
    i = h.find("(")
    fixed = h[:i] if i >= 0 else h
    total, atoms = 0, set()
    for cnt, _sign, alist in MOBILE.findall(h[i:] if i >= 0 else ""):
        total += int(cnt) if cnt else 1
        atoms.update(int(a) for a in alist.split(",") if a)
    return fixed.rstrip(","), total, frozenset(atoms)


def main():
    cz = sqlite3.connect(str(require(ZENODO_CACHE, "upstream Zenodo cache")))
    zrows = defaultdict(list)
    for cid, ik, inchi, dc in cz.execute(
            "select id,inchi_key,inchi,dissociation_constants from compounds "
            "where inchi is not null"):
        L = split_inchi(inchi)
        if not L:
            continue
        v = pickle.loads(dc) if isinstance(dc, (bytes, bytearray)) else dc
        zrows[(L["formula"], L["c"])].append((cid, ik, L, tuple(v or [])))
    zkeys = {ik for rows in zrows.values() for _, ik, _, _ in rows if ik}

    ours = {}
    with require(STRUCTURES, "ModelSEED structures").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Type"] in ("InChI", "InChIKey"):
                ours.setdefault(r["ID"], {})[r["Type"]] = r["Structure"]

    rank = list(csv.DictReader(
        require(RANKED, "polyprotic ranking").open(), delimiter="\t"))

    verdicts = defaultdict(lambda: [0, 0])
    multi = 0
    detail = defaultdict(list)
    pdiff = 0
    for r in rank:
        sid, n = r["seed_id"], int(r["n_reactions"])
        o = ours.get(sid, {})
        key, L = o.get("InChIKey"), split_inchi(o.get("InChI"))
        if not L or (key and key in zkeys):
            continue                       # absent, or already an exact match
        # Only candidates that agree on stereo; the stereo-differing class is
        # handled separately.
        cands = [c for c in zrows.get((L["formula"], L["c"]), [])
                 if c[3] and c[2]["t"] == L["t"]]
        if not cands:
            continue
        if len(cands) > 1:
            multi += 1

        # Selection must SEARCH the candidates, not take the first. IPP and
        # DMAPP are constitutional isomers differing only in /h, and Zenodo
        # holds a row for each; taking candidates[0] compares DMAPP against the
        # IPP row and rejects a mismatch that selection itself manufactured.
        of, om, oa = parse_h(L["h"])
        exact = [c for c in cands if c[2]["h"] == L["h"]]
        regrouped = [c for c in cands if parse_h(c[2]["h"]) == (of, om, oa)]
        if exact:
            v = "safe: /h identical (differs only in /p)"
            cid, zik, ZL, lad = exact[0]
        elif regrouped:
            v = "safe: same protons, different grouping"
            cid, zik, ZL, lad = regrouped[0]
        else:
            cid, zik, ZL, lad = cands[0]
            zf, zm, za = parse_h(ZL["h"])
            v = ("UNSAFE: fixed-H layer differs" if of != zf else
                 "UNSAFE: mobile proton COUNT differs" if om != zm else
                 "UNSAFE: mobile proton ATOM SET differs")
        if L["p"] != ZL["p"]:
            pdiff += 1
        verdicts[v][0] += 1
        verdicts[v][1] += n
        detail[v].append((n, sid, r["name"][:30], L["p"] or "p0", ZL["p"] or "p0"))

    tot = sum(v[1] for v in verdicts.values())
    print("CONNECTIVITY+STEREO MATCHES WHERE THE InChIKey DIFFERS\n")
    for k, (c, rx) in sorted(verdicts.items(), key=lambda x: -x[1][1]):
        print(f"  {k:44s} {c:4,} cpds {rx:6,} rxns ({100*rx/tot:4.1f}%)")
    print(f"\n  of these, /p (protonation offset) differs in {pdiff:,}")
    print(f"  compounds with >1 candidate on the same formula+/c+/t: {multi:,}"
          f"  (selection must search these, not take the first -- see IPP/DMAPP)")
    for k in sorted(detail):
        if k.startswith("UNSAFE"):
            print(f"\n-- {k}:")
            for n, s, nm, op, zp in sorted(detail[k], reverse=True)[:12]:
                print(f"     {s} {n:>4} {nm:32s} ours/{op} zenodo/{zp}")
    print("\n-- sample of the safe-regrouping class:")
    for n, s, nm, op, zp in sorted(detail["safe: same protons, different grouping"],
                                   reverse=True)[:8]:
        print(f"     {s} {n:>4} {nm:32s} ours/{op} zenodo/{zp}")

    # Emit the classification so the evidence grader can consume it rather than
    # re-deriving it, and so the refused set is a reviewable artefact.
    DATA_OUT.mkdir(parents=True, exist_ok=True)
    out = DATA_OUT / "structural_match_classification.tsv"
    with out.open("w") as fh:
        fh.write("seed_id\tname\tn_reactions\tverdict\tsafe\tours_p\tzenodo_p\n")
        for k, entries in sorted(detail.items()):
            safe = 0 if k.startswith("UNSAFE") else 1
            for n, s, nm, op, zp in sorted(entries, reverse=True):
                fh.write(f"{s}\t{nm}\t{n}\t{k}\t{safe}\t{op}\t{zp}\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
