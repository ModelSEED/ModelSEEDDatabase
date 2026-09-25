#!/usr/bin/env python
"""Rebuild the reaction->pathway annotation for every reaction, old and new.

WHY THIS EXISTS. Biochemistry/Aliases/Unique_ModelSEED_Reaction_Pathways.txt
stopped at rxn48568 and the per-reaction `pathways` field stopped with it, so
not one of the 12,261 reactions added since the 2020 release carried a pathway
annotation -- in the alias file or in the record. The step that used to produce
it was Scripts/Archived_Perl_Scripts/Compile_External_Pathways.pl, which was
never ported when the pipeline moved to Python.

THIS IS ADDITIVE, NOT A REGENERATION, and the distinction matters. The
committed pathway tables are a smaller snapshot than the MetaCyc and KEGG
distributions that produced the original alias file: rebuilding every row from
them reproduces only 104,203 of the 121,444 rows on disk and would silently drop
17,241 -- every row of 601 reactions among them. Those distributions are
licence-restricted and are not in this repository, so a faithful regeneration
is not possible from a clean checkout and is NOT attempted here.

What this does instead is leave every existing row untouched and annotate only
reactions that currently carry none. That is a strict improvement with no
regression. It does not close the gap -- the reactions it cannot reach still
need the real pipeline port against the upstream sources -- and --check
reports exactly how many remain.

WHAT IT READS (all committed; no downloads):
    Scripts/Provenance/MetaCyc/MetaCyc_pathways.tsv
    Scripts/Provenance/KEGG/KEGG_pathways.tsv
    Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt

WHAT IT WRITES:
    Biochemistry/Aliases/Unique_ModelSEED_Reaction_Pathways.txt
    the `pathways` field of Biochemistry/reaction_NN.{json,tsv}

BOTH pathway tables mix individual pathways with ontology CLASSES in one id
column -- `Degradation`, `Biosynthesis` and `Energy-Metabolism` are ids, not
just parents. A reaction is annotated with its own pathways AND the transitive
closure of their parents, which is what the pre-2020 rows contain and what
makes a class-level distribution possible.

    ~/Documents/py_venv/bin/python Build_Reaction_Pathways.py [--check]

--check compares what would be written against what is on disk and writes
nothing. Use it to confirm that a rebuild preserves every pre-existing row
before letting it touch the release.
"""
import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ALIASES = ROOT / "Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt"
OUT = ROOT / "Biochemistry/Aliases/Unique_ModelSEED_Reaction_Pathways.txt"
TABLES = {
    "MetaCyc": ROOT / "Scripts/Provenance/MetaCyc/MetaCyc_pathways.tsv",
    "KEGG": ROOT / "Scripts/Provenance/KEGG/KEGG_pathways.tsv",
}

# Which ids each source is allowed to emit. KEGG's table is a module table:
# its rows are M-numbers and its parent column mixes rn pathway ids with
# free-text category names. The shipped annotation contains rn ids only, so
# anything else is dropped rather than emitted with an empty name.
KEEP = {"KEGG": lambda i: i.startswith("rn"), "MetaCyc": lambda i: True}
SOURCE_ORDER = ("MetaCyc", "KEGG")


def load_table(path):
    """id -> (name, parents, reactions). Classes and pathways share the column."""
    name, parents, rxns = {}, {}, defaultdict(set)
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ident = r["id"]
            name[ident] = (r.get("name") or "").strip()
            parents[ident] = [p for p in (r.get("parent") or "").split("|") if p]
            for x in (r.get("reactions") or "").split("|"):
                if x:
                    rxns[x].add(ident)
    return name, parents, rxns


def closure(idents, parents, limit=32):
    """Every ancestor of every id, transitively. `limit` guards a cyclic
    ontology -- MetaCyc's is a DAG in principle but this must not hang a
    release build on a bad edit upstream."""
    seen, frontier, depth = set(idents), set(idents), 0
    while frontier and depth < limit:
        nxt = set()
        for i in frontier:
            for p in parents.get(i, ()):
                if p not in seen:
                    seen.add(p)
                    nxt.add(p)
        frontier, depth = nxt, depth + 1
    return seen


def seed_to_external():
    """ModelSEED reaction id -> {source: {external ids}}."""
    m = defaultdict(lambda: defaultdict(set))
    with ALIASES.open() as fh:
        rd = csv.reader(fh, delimiter="\t"); next(rd)
        for row in rd:
            if len(row) >= 3 and row[2] in TABLES:
                m[row[0]][row[2]].add(row[1])
    return m


def build():
    """rxn -> source -> sorted ['ID (Name)'] strings."""
    tables = {src: load_table(p) for src, p in TABLES.items()}
    known = known_names()
    per = defaultdict(lambda: defaultdict(set))
    for seed, by_src in seed_to_external().items():
        for src, externals in by_src.items():
            name, parents, rxns = tables[src]
            direct = {p for ext in externals for p in rxns.get(ext, ())}
            if not direct:
                continue
            keep = KEEP[src]
            for ident in closure(direct, parents):
                if not keep(ident) or ident not in name:
                    continue
                label = name[ident] or known.get(ident, "")
                per[seed][src].add(f"{ident} ({label})")
    return per


def known_names():
    """Display names already in the alias file, harvested as "ID (Name)".

    683 of the 4,128 rows in MetaCyc_pathways.tsv are ontology classes with an
    EMPTY name column -- Degradation, Antibiotic-Biosynthesis, Fermentation and
    the rest. The pre-2020 rows carry proper names for most of them because
    they were built from the full MetaCyc distribution, which is not in this
    repository. Rather than invent a name or emit "Antibiotic-Biosynthesis ()"
    12,913 times, reuse what the file already knows. Where nothing is known the
    empty form is kept, which is the convention already present on 937 rows.
    """
    names = {}
    with OUT.open() as fh:
        rd = csv.reader(fh, delimiter="\t"); next(rd)
        for row in rd:
            if len(row) < 2 or " (" not in row[1] or not row[1].endswith(")"):
                continue
            ident, _, label = row[1].partition(" (")
            label = label[:-1]
            if label and ident not in names:
                names[ident] = label
    return names


def existing():
    """Every row currently in the alias file, and the reactions it covers."""
    rows, covered = set(), set()
    with OUT.open() as fh:
        rd = csv.reader(fh, delimiter="\t"); next(rd)
        for row in rd:
            if len(row) >= 3:
                rows.add((row[0], row[1], row[2]))
                covered.add(row[0])
    return rows, covered


def additive(per):
    """Existing rows, plus rebuilt rows for reactions that have none."""
    rows, covered = existing()
    added = {(s, e, src) for s, by in per.items() if s not in covered
             for src, es in by.items() for e in es}
    return rows | added, added


def write_alias_file(rows):
    with OUT.open("w") as fh:
        fh.write("ModelSEED ID\tExternal ID\tSource\n")
        for seed, ext, src in sorted(rows):
            fh.write(f"{seed}\t{ext}\t{src}\n")
    return len(rows)


def write_records(rows):
    """Set the `pathways` field for records that have none, matching the
    shipped format exactly: one string per source, with the ids sorted --

        "MetaCyc: Degradation (Degradation/Utilization/Assimilation); ..."

    Records that already carry a value are left exactly as they are, for the
    same reason the alias file is not regenerated."""
    by_rxn = defaultdict(lambda: defaultdict(set))
    for seed, ext, src in rows:
        by_rxn[seed][src].add(ext)
    touched = 0
    for shard in sorted(ROOT.glob("Biochemistry/reaction_[0-9][0-9].json")):
        data = json.load(open(shard))
        changed = False
        for rxn in data:
            if rxn.get("pathways"):
                continue
            by_src = by_rxn.get(rxn["id"])
            if not by_src:
                continue
            # Source order matches every existing multi-source record, all
            # 2,023 of which are MetaCyc first, then KEGG -- not alphabetical.
            rxn["pathways"] = [f"{src}: " + "; ".join(sorted(by_src[src]))
                               for src in SOURCE_ORDER if src in by_src]
            changed = True
            touched += 1
        if changed:
            shard.write_text(json.dumps(data, indent=4, sort_keys=True))
            _sync_tsv(shard.with_suffix(".tsv"), data)
    return touched


def _sync_tsv(path, data):
    """Mirror the pathways column into the flat file, leaving every other byte
    alone.

    Deliberately NOT csv.writer. These files are raw tab-delimited text, not
    RFC-4180: the stoichiometry column contains unquoted " characters, so a
    csv round-trip quotes every such field and rewrites all 50 files for a
    change that touches a few hundred lines."""
    by_id = {r["id"]: r for r in data}
    lines = path.read_text().split("\n")
    if not lines:
        return
    header = lines[0].split("\t")
    if "pathways" not in header:
        return
    col = header.index("pathways")
    for i, line in enumerate(lines[1:], start=1):
        if not line:
            continue
        cells = line.split("\t")
        if len(cells) <= col:
            continue
        rxn = by_id.get(cells[0])
        if rxn is None:
            continue
        pw = rxn.get("pathways")
        cells[col] = "|".join(pw) if pw else "null"
        lines[i] = "\t".join(cells)
    path.write_text("\n".join(lines))


def check(per):
    """Report what an additive run would change, and what it leaves undone."""
    have, covered = existing()
    merged, added = additive(per)
    live = set()
    for shard in sorted(ROOT.glob("Biochemistry/reaction_[0-9][0-9].json")):
        live |= {r["id"] for r in json.load(open(shard))}
    gained = {r for r, _, _ in added}
    print(f"alias rows on disk        {len(have):>9,}")
    print(f"  preserved               {len(have & merged):>9,}   (must equal the line above)")
    print(f"  added                   {len(added):>9,}")
    print(f"reactions annotated       {len(covered):>9,}  ->  {len(covered | gained):,}")
    print(f"  newly annotated         {len(gained):>9,}")
    print(f"still unannotated         {len(live - covered - gained):>9,}"
          f"   of {len(live):,} reactions")
    full = {s for s, by in per.items() for _ in by}
    print(f"\nFor reference, a full rebuild from the committed tables would")
    print(f"reproduce only {len(have & {(s, e, sr) for s, by in per.items() for sr, es in by.items() for e in es}):,}"
          f" of those {len(have):,} rows, which is why this is additive.")
    return len(have & merged) == len(have)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="compare against what is on disk; write nothing")
    args = ap.parse_args()
    per = build()
    if args.check:
        return 0 if check(per) else 1
    merged, added = additive(per)
    if not added:
        print("nothing to add; alias file already covers every reaction it can")
        return 0
    n = write_alias_file(merged)
    t = write_records(merged)
    print(f"wrote {n:,} rows to {OUT.relative_to(ROOT)}")
    print(f"updated the pathways field on {t:,} reaction records")
    return 0


if __name__ == "__main__":
    sys.exit(main())
