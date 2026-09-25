#!/usr/bin/env python
"""Where the 2020->2026 growth landed, at MetaCyc class ("subsystem") level.

Written for reviewer 1 comment 1: "it would be curious to see where are we
still gaining new information". The answer is a frequency distribution over
MetaCyc pathway classes, computed for the reactions that are new since the
2020 release against the ones that were already there, so the two can be read
side by side.

Baseline is commit fd6c7849 (2020-11-10), the last commit of 2020 and the same
baseline figure_common._growth() uses for the compound panel -- do not swap it
for a different 2020 commit without changing that function too.

POPULATION. Counts every record, obsolete included, because that is the basis
the manuscript uses everywhere else -- the abstract's "~56,000 reactions", Fig
2A's 56,002 and M12's 25,855 are all all-records figures. On that basis growth
is +28.0%, NOT the +55% M09 claims: 55% only appears if 2026 totals are set
against a 2020 baseline with its 7,577 obsolete rows removed. Pass --live for
the live-vs-live basis (+33.7%), which is defensible but would require
converting every other number in the paper with it.

Run:  ~/Documents/py_venv/bin/python pathway_distribution_of_growth.py
"""
import csv
import io
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
C2020 = "fd6c7849891ef4bbeb6eac072f5a6f7adff05b0e"
RXN_2020 = "Biochemistry/reactions.tsv"
PATHWAYS = ROOT / "Biochemistry/Aliases/Unique_ModelSEED_Reaction_Pathways.txt"
ECS = ROOT / "Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt"
MC_PWY = ROOT / "Scripts/Provenance/MetaCyc/MetaCyc_pathways.tsv"
RXN_SRC = ROOT / "Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt"


def _rows(text):
    return csv.DictReader(io.StringIO(text), delimiter="\t")


LIVE_ONLY = "--live" in sys.argv


def live_2020():
    """Reaction IDs in the 2020 release (all records unless --live)."""
    out = subprocess.run(["git", "show", f"{C2020}:{RXN_2020}"], cwd=str(ROOT),
                         capture_output=True, text=True, check=True).stdout
    return {r["id"] for r in _rows(out)
            if not (LIVE_ONLY and r.get("is_obsolete") == "1")}


def live_2026():
    """Non-obsolete reaction IDs in the current release, across the shards."""
    ids = set()
    for shard in sorted(ROOT.glob("Biochemistry/reaction_[0-9][0-9].tsv")):
        with open(shard) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if not (LIVE_ONLY and r.get("is_obsolete") == "1"):
                    ids.add(r["id"])
    return ids


def metacyc_tables():
    """Parse MetaCyc_pathways.tsv into (reaction -> pathways, parent map, names,
    class ids).

    CAREFUL: that file is not a list of pathways. It mixes individual pathways
    (PWY-7805) with ontology CLASSES (Degradation, Biosynthesis,
    Energy-Metabolism), all in the same id column -- an earlier version of this
    script assumed every id was a pathway and consequently found zero classes.
    A class is identified as any id that some other row names as its parent.
    """
    rx2pwy, parent, names = defaultdict(set), {}, {}
    with open(MC_PWY) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            names[r["id"]] = r["name"] or r["id"]
            parent[r["id"]] = [p for p in (r.get("parent") or "").split("|") if p]
            for rx in (r["reactions"] or "").split("|"):
                if rx:
                    rx2pwy[rx].add(r["id"])
    classes = {p for ps in parent.values() for p in ps}
    return rx2pwy, parent, names, classes


def seed_to_metacyc():
    """ModelSEED reaction id -> MetaCyc reaction ids. This alias file DOES run
    to rxn60859, which is why the post-2020 intake can be re-annotated at all."""
    m = defaultdict(set)
    with open(RXN_SRC) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Source"] == "MetaCyc":
                m[r["ModelSEED ID"]].add(r["External ID"])
    return m


def classes_for(seed_ids, s2m, rx2pwy, parent, classes):
    """Distinct-reaction counts per MetaCyc class, walking each pathway up to
    its parents. Counts REACTIONS, not (reaction, pathway) pairs -- a reaction
    in three pathways of one class must count once."""
    per = Counter()
    reached = set()
    for sid in seed_ids:
        cls = set()
        for mc in s2m.get(sid, ()):
            for pwy in rx2pwy.get(mc, ()):
                cls |= {p for p in parent.get(pwy, []) if p in classes or True}
        if cls:
            reached.add(sid)
        for c in cls:
            per[c] += 1
    return per, reached


def reaction_to_classes(pathway_ids):
    """rxn id -> set of MetaCyc class ids (ancestors), and -> set of pathway ids."""
    classes, pathways = defaultdict(set), defaultdict(set)
    names = {}
    with open(PATHWAYS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Source"] != "MetaCyc":
                continue
            ext = r["External ID"]
            ident, _, label = ext.partition(" (")
            label = label.rstrip(")") or ident
            names[ident] = label
            if ident in pathway_ids:
                pathways[r["ModelSEED ID"]].add(ident)
            else:
                classes[r["ModelSEED ID"]].add(ident)
    return classes, pathways, names


def ec_class(rxn_ids):
    """First digit of the EC number -> broad enzyme class, as a fallback for
    reactions MetaCyc never placed in a pathway (most of the Rhea intake)."""
    EC = {"1": "oxidoreductase", "2": "transferase", "3": "hydrolase",
          "4": "lyase", "5": "isomerase", "6": "ligase", "7": "translocase"}
    have = Counter()
    seen = set()
    with open(ECS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rid = r["ModelSEED ID"]
            if rid not in rxn_ids or rid in seen:
                continue
            top = r["External ID"].split(".")[0]
            if top in EC:
                have[EC[top]] += 1
                seen.add(rid)
    return have, len(seen)


def source_of(rxn_ids):
    """Which primary database supplies each reaction (a reaction may have more
    than one, so these do not sum to the set size)."""
    per = defaultdict(set)
    with open(RXN_SRC) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["ModelSEED ID"] in rxn_ids:
                per[r["Source"]].add(r["ModelSEED ID"])
    return per


def table(title, new_c, old_c, n_new, n_old, names, limit=22):
    print(f"\n{title}")
    print(f"{'class':<52}{'new':>7}{'% new':>8}{'old':>8}{'% old':>8}{'ratio':>8}")
    print("-" * 91)
    for ident, n in new_c.most_common(limit):
        o = old_c.get(ident, 0)
        pn, po = 100 * n / n_new, 100 * o / n_old if n_old else 0
        ratio = f"{pn / po:5.2f}" if po else "  inf"
        print(f"{names.get(ident, ident)[:51]:<52}{n:>7}{pn:>7.1f}%{o:>8}{po:>7.1f}%{ratio:>8}")


def main():
    old, new_all = live_2020(), live_2026()
    gained = new_all - old
    kept = new_all & old
    print("=" * 91)
    print(f"REACTION SETS ({'live only (--live)' if LIVE_ONLY else 'ALL records -- the manuscript basis'})")
    print("=" * 91)
    print(f"2020 reactions               {len(old):>8,}")
    print(f"2026 reactions               {len(new_all):>8,}")
    print(f"  of which new since 2020    {len(gained):>8,}")
    print(f"  of which carried over      {len(kept):>8,}")
    print(f"  2020 ids now absent        {len(old - new_all):>8,}")
    print(f"growth on this basis         {100 * (len(new_all) / len(old) - 1):>7.1f}%")

    print("\nprimary source of the new reactions (sets overlap):")
    for src, ids in sorted(source_of(gained).items(), key=lambda kv: -len(kv[1])):
        if len(ids) >= 200:
            print(f"  {src:<24}{len(ids):>8,}")

    rx2pwy, parent, names, classes = metacyc_tables()
    s2m = seed_to_metacyc()

    # Both sides go through the SAME join, so the two columns are comparable.
    # Reading the old side out of Unique_ModelSEED_Reaction_Pathways.txt instead
    # would compare a full ancestor closure against a direct-parent lookup.
    new_c, new_reached = classes_for(gained, s2m, rx2pwy, parent, classes)
    old_c, old_reached = classes_for(kept, s2m, rx2pwy, parent, classes)
    print(f"\nreachable through the MetaCyc pathway join: "
          f"{len(new_reached):,} of {len(gained):,} new, "
          f"{len(old_reached):,} of {len(kept):,} carried over")
    table("MetaCyc CLASS distribution, new vs carried-over, DISTINCT REACTIONS "
          "(ratio > 1 = enriched among the new)",
          new_c, old_c, len(new_reached), len(old_reached), names)

    ENERGY = {"Energy-Metabolism", "Respiration", "Fermentation", "Photosynthesis",
              "Electron-Transfer", "TCA-VARIANTS", "Glycolysis",
              "Pentose-Phosphate-Cycle", "Chemoautotrophic-Energy-Metabolism",
              "Methanogenesis"}
    def _in_energy(ids):
        n = 0
        for sid in ids:
            cls = {p for mc in s2m.get(sid, ())
                   for pwy in rx2pwy.get(mc, ()) for p in parent.get(pwy, [])}
            n += bool(cls & ENERGY)
        return n
    print(f"\ncentral/energy metabolism classes {sorted(ENERGY)[:3]}...:")
    print(f"  new reactions landing there   {_in_energy(gained):>7,} "
          f"of {len(gained):,} new")
    print(f"  carried-over reactions there  {_in_energy(kept):>7,} "
          f"of {len(kept):,}")

    classes_by_rxn, pathways, _n2 = reaction_to_classes(
        {r["id"] for r in csv.DictReader(open(MC_PWY), delimiter="\t")})
    ann_new = sum(1 for r in gained if classes_by_rxn.get(r) or pathways.get(r))
    ann_old = sum(1 for r in kept if classes_by_rxn.get(r) or pathways.get(r))
    print(f"\nMetaCyc pathway annotation coverage")
    print(f"  new reactions with any MetaCyc pathway/class "
          f"{ann_new:>7,} / {len(gained):,}  ({100 * ann_new / len(gained):.1f}%)")
    print(f"  carried-over reactions with any              "
          f"{ann_old:>7,} / {len(kept):,}  ({100 * ann_old / len(kept):.1f}%)")

    unann = {r for r in gained if not (classes_by_rxn.get(r) or pathways.get(r))}
    ec_new, n_ec = ec_class(unann)
    print(f"\nEC fallback for the {len(unann):,} new reactions MetaCyc never placed:")
    print(f"  carry an EC number   {n_ec:>7,}  ({100 * n_ec / len(unann):.1f}%)")
    for k, v in ec_new.most_common():
        print(f"    {k:<20}{v:>7,}  ({100 * v / n_ec:.1f}% of EC-annotated)")

    # Via the reconstruction, NOT the alias file -- the alias file holds no
    # post-2020 reaction at all, so reading it here would silently print
    # nothing and read as "no pathways gained" rather than "never annotated".
    most = Counter()
    for sid in gained:
        for mc in s2m.get(sid, ()):
            for pwy in rx2pwy.get(mc, ()):
                if pwy not in classes:      # individual pathways only
                    most[pwy] += 1
    print(f"\nMost-gained individual MetaCyc pathways (new reactions only, "
          f"via the reconstruction):")
    for ident, n in most.most_common(15):
        print(f"  {names.get(ident, ident)[:60]:<62}{n:>6}")


if __name__ == "__main__":
    sys.exit(main())
