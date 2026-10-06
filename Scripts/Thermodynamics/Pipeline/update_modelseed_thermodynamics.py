"""Write the regenerated eQuilibrator energies into ModelSEED's JSON.

Touches exactly one key -- thermodynamics["eQuilibrator"] -- in
Biochemistry/{reaction,compound}_*.json. Everything else is left alone:
top-level deltag/deltagerr (those are the canonical promoted values, and
promotion is a separate step that must run after reversibility), the .tsv
files, and every other thermodynamics source.

Stored format is [deltag, deltagerr, reversibility] in **kcal/mol**.

Reversibility is NOT recomputed. An existing flag is carried through
untouched; an entry that is newly covered has no prior flag and gets "?",
the placeholder already used 28,689 times in this database. Those ids are
written to a report so the reversibility pass knows exactly what to look at.

An entry that had an eQuilibrator value and has none now loses the key. The
old values are preserved in eQuilibrator/data/eQuilibrator-2020/ -- leaving a
stale entry beside fresh ones under the same name would be a trap.

Round-trips byte-identically with indent=4, so the diff shows only real
changes. Pass --dry-run to see the counts without writing.
"""

import os
import argparse
import csv
import re
import hashlib
import glob
import json
import sys
from collections import Counter
from pathlib import Path

# Relocated from the eQuilibrator working tree 2026-09-05. ROOT is that
# tree -- it holds the caches and fitted parameters, which are gigabytes and
# are not in this repository -- so it is named by environment variable rather
# than derived from this file's location.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))
# Default target: the main ModelSEEDDatabase clone. Override with --worktree to
# write into a git worktree or an integration branch checkout instead. This used
# to be pinned to the eq-rerun worktree, which meant the script silently wrote
# to the wrong tree once that branch was superseded.
# This repository, derived from the script location now that it lives in-tree.
DEFAULT_MSD = Path(__file__).resolve().parents[3]
SOURCE = "eQuilibrator"
UNKNOWN_REVERSIBILITY = "?"


PROVENANCE = re.compile(r"(\w+)=(\S+)")
CONDITIONS = ("p_h", "ionic_strength", "p_mg", "T")


def read_provenance(path):
    """Parse the '# cache=... params=... p_h=...' line the generators write."""
    with open(path) as fh:
        first = fh.readline()
    if not first.startswith("#") or "cache=" not in first:
        raise SystemExit(
            f"{path}: no provenance header. Every generated energy table starts "
            f"'# cache=... params=... p_h=...'; a file without one did not come "
            f"from generate_modelseed_*_energies.py and must not be ingested.")
    return dict(PROVENANCE.findall(first))


def cache_fingerprint(declared):
    """Content hash of the cache a table declares, or None if it is gone.

    Compared by CONTENT, not by the recorded path: a cache is routinely
    promoted under a new name (cache_win -> cache_final) between generating a
    table and shipping it, so path equality would reject a legitimate run,
    while content equality catches the case that matters -- a table generated
    from a different cache than the one being published.
    """
    path = Path(declared)
    if not path.is_absolute():
        path = ROOT / path
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_inputs(paths, shipped_cache, allow_stale=False):
    """Refuse to ingest tables that did not come from the cache being shipped.

    The failure this exists to prevent is silent and looks healthy: on
    2026-09-04 both input tables were five days stale, the run reported "23,800
    updated" and the entry count rose, and the JSON received pre-rebuild
    numbers. Nothing in the summary indicated it. A wrong ingest must be loud.
    """
    provs = {p: read_provenance(p) for p in paths}

    for key in CONDITIONS:
        seen = {prov.get(key) for prov in provs.values()}
        if len(seen) > 1:
            raise SystemExit(
                f"input tables disagree on {key}: {sorted(map(str, seen))}. "
                f"They were generated under different conditions and cannot be "
                f"ingested together.")

    fps = {}
    for path, prov in provs.items():
        fp = cache_fingerprint(prov["cache"])
        if fp is None:
            msg = (f"{path} names cache {prov['cache']}, which no longer exists, "
                   f"so its provenance cannot be verified.")
            if not allow_stale:
                raise SystemExit(msg + " Re-generate it, or pass --allow-stale.")
            print(f"WARNING: {msg}")
        fps[path] = fp

    shipped = cache_fingerprint(shipped_cache)
    if shipped is None:
        raise SystemExit(f"shipped cache not found: {shipped_cache}")

    bad = [f"{p} (from {provs[p]['cache']})" for p, fp in fps.items()
           if fp is not None and fp != shipped]
    if bad:
        msg = ("input tables were NOT generated from the cache being shipped "
               f"({shipped_cache}):\n  " + "\n  ".join(bad))
        if not allow_stale:
            raise SystemExit(
                msg + "\n\nRe-run the generators against that cache, or pass "
                "--allow-stale if you genuinely mean to publish tables from a "
                "different one.")
        print(f"WARNING: {msg}")

    cond = next(iter(provs.values()))
    print("inputs verified against " + str(shipped_cache))
    print("  conditions: " + "  ".join(f"{k}={cond.get(k)}" for k in CONDITIONS))


def load_new(path, id_col, dg_col, err_col):
    fh = open(path)
    fh.readline()                      # '# cache=... p_mg=...' provenance line
    out = {}
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["status"] != "ok":
            continue
        try:
            out[r[id_col]] = (round(float(r[dg_col]), 3), round(float(r[err_col]), 3))
        except (ValueError, KeyError):
            pass
    return out


def resolve_biochemistry(root):
    """Accept either a repo root or its Biochemistry/ directory.

    Verified rather than assumed: a mistyped path would otherwise glob to
    nothing and the script would report a clean run having written nothing.
    """
    root = Path(root).expanduser().resolve()
    bio = root if root.name == "Biochemistry" else root / "Biochemistry"
    if not (bio / "reaction_00.json").is_file():
        sys.exit(f"not a ModelSEEDDatabase checkout: {root}\n"
                 f"  (expected {bio / 'reaction_00.json'})")
    return bio


def apply(kind, new, dry, wt):
    # Reactions store [dg, err, reversibility]; compounds store [dg, err].
    # A formation energy has no direction, so writing a reversibility slot
    # into a compound would invent a field the schema does not have.
    triple = kind == "reaction"
    stats = Counter()
    newly_covered = []
    for f in sorted(glob.glob(str(wt / f"{kind}_*.json"))):
        data = json.load(open(f, encoding="utf-8"))
        changed = False
        for rec in data:
            thermo = rec.get("thermodynamics")
            rid = rec["id"]
            fresh = new.get(rid)
            if fresh is None:
                if thermo and SOURCE in thermo:
                    del thermo[SOURCE]
                    stats["removed (no new value)"] += 1
                    changed = True
                else:
                    stats["untouched (no value either way)"] += 1
                continue
            if thermo is None:
                thermo = rec["thermodynamics"] = {}
            prior = thermo.get(SOURCE)
            if not triple:
                thermo[SOURCE] = [fresh[0], fresh[1]]
                stats["updated" if prior else "added"] += 1
                changed = True
                continue
            if prior and len(prior) > 2:
                rev = prior[2]
                stats["updated (reversibility carried over)"] += 1
            else:
                rev = UNKNOWN_REVERSIBILITY
                stats["added (reversibility unknown)"] += 1
                newly_covered.append(rid)
            thermo[SOURCE] = [fresh[0], fresh[1], rev]
            changed = True
        if changed and not dry:
            with open(f, "w", encoding="utf-8") as fh:
                fh.write(json.dumps(data, indent=4) + "\n")
    return stats, newly_covered


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--cache", type=Path,
                    default=ROOT / "data/cache_final/compounds.sqlite",
                    help="the cache being shipped. Input tables must declare a "
                         "cache with identical content. Default: %(default)s")
    ap.add_argument("--allow-stale", action="store_true",
                    help="downgrade the provenance mismatch to a warning. Only "
                         "for deliberately publishing tables from another cache.")
    ap.add_argument("--worktree", type=Path, default=DEFAULT_MSD,
                    help="ModelSEEDDatabase checkout to write into (repo root "
                         "or its Biochemistry/ dir). Default: %(default)s")
    args = ap.parse_args()

    wt = resolve_biochemistry(args.worktree)
    print(f"target: {wt}")

    rxn_tsv = ROOT / "data/modelseed_energies.tsv"
    cpd_tsv = ROOT / "data/modelseed_formation_energies.tsv"
    check_inputs([rxn_tsv, cpd_tsv], args.cache, allow_stale=args.allow_stale)

    rxn = load_new(rxn_tsv, "reaction_id",
                   "dg_prime_kcal_per_mol", "uncertainty_kcal_per_mol")
    cpd = load_new(cpd_tsv, "compound_id",
                   "dgf_prime_kcal_per_mol", "uncertainty_kcal_per_mol")
    print(f"new values: {len(rxn)} reactions, {len(cpd)} compounds"
          f"{'   [DRY RUN]' if args.dry_run else ''}\n")

    report = []
    for kind, new in (("reaction", rxn), ("compound", cpd)):
        stats, fresh_ids = apply(kind, new, args.dry_run, wt)
        print(f"{kind}:")
        for k, v in stats.most_common():
            print(f"    {k:<40}{v:7d}")
        report += [(kind, i) for i in fresh_ids]

    if report and not args.dry_run:
        p = ROOT / "data/newly_covered_needs_reversibility.tsv"
        with p.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["kind", "id"])
            w.writerows(report)
        print(f"\nwrote {p}  ({len(report)} entries need a reversibility call)")


if __name__ == "__main__":
    main()
