"""Build a compound cache that contains ModelSEED and nothing else.

The shipped Zenodo cache is 699,184 compounds -- MetaNetX's universe, half of
it SwissLipids -- of which only ~36,800 ever carry a seed: accession and ~640
affect the fit. `build_path_a_cache.py` overlays our accessions onto that,
which works but leaves a `retain` fallback: where we cannot supply a structure
it puts the MetaNetX-derived mapping back, so some seed: accessions still
resolve to somebody else's molecule.

This builds the cache the other way round. Start empty; put in only what
ModelSEED needs. With no prior mapping there is nothing to retain, so the
property "every seed: accession resolves to a structure we curate" is
structural rather than a policy setting.

Three populations go in:

  carry-over    22,326 compounds whose structure the Zenodo cache already
                holds (`our_cache_id` set). Copying that row is NOT a
                provenance compromise -- the cached compound *is* our
                molecule, which is what the match means. This path is
                load-bearing, not an optimisation: phosphate, water, CO2 and
                ammonia cannot be created at all, because group decomposition
                rejects small inorganics, and several of them are anchors.

  created       6,213 compounds whose structure the cache lacks, built from
                our SMILES + pKas.

  training-only the accessions the training data references that have no
                ModelSEED route -- metanetx.chemical: and coco:, plus the
                kegg: that do not resolve. Unavoidable: train_path_b resolves
                training reactions through ccache.get_compound, so these must
                exist or the training set silently shrinks. Copied verbatim
                and reported separately, because they are the one place a
                foreign structure legitimately remains: the measurements are
                keyed to KEGG and that cannot be wished away.

Not represented, and reported rather than papered over: 1,839 compounds whose
structure we hold but that the cache lacks and we have no pKas to build
(fixable -- 6,405 ModelSEED compounds have simply never been through MolGpKa),
and 6,565 where we hold no structure at all.

Defect carried forward: `pick_compound` (prefer an anchor, else lowest id).
1,856 InChIKeys in the source cache map to more than one compound, and naively
taking the highest id demoted water, phosphate, CO2, ammonia and pyruvate onto
unusable duplicates -- passing every integrity check while making every
reaction touching them infinitely uncertain.

Defect replaced: the old builder's `retain` fallback existed so a failed
creation would not silently drop an accession. There is no prior mapping to
fall back to here, so failures are written to a report and the run exits
non-zero past --max-failures. Loud beats plausible.
"""

import argparse
import csv
import os
import pickle
import sqlite3
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

# Relocated from the eQuilibrator working tree 2026-09-05. ROOT is that
# tree -- it holds the caches and fitted parameters, which are gigabytes and
# are not in this repository -- so it is named by environment variable rather
# than derived from this file's location.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

SOURCE_CACHE = Path(os.path.expanduser("~/.cache/equilibrator/compounds.sqlite"))
MAPPING = ROOT / "data" / "seed_mapping.tsv"
DEFAULT_OUT = ROOT / "data" / "modelseed_cache_msd" / "compounds.sqlite"
FAILURES = ROOT / "data" / "cache_creation_failures.tsv"
SEED_NAMESPACE = "seed"


def classify(row, must_carry=None, predicted=None):
    """carry | create | blocked | no_structure, plus the reason.

    ``must_carry`` switches on no-carry-over mode: a compound whose structure
    the Zenodo cache happens to hold is still BUILT from our own structure
    unless it is in that set. Carrying a Zenodo row over is not provenance
    neutral -- equilibrator-assets builds that cache by shelling out to
    ChemAxon cxcalc, so every carried row imports a Marvin-derived microspecies
    even when the molecule is ours.

    The set is not small and is not optional. eQuilibrator's group decomposer
    refuses 3,526 of the 22,326 carry-overs outright (measured by
    tools/must_carry_over.py), and the handful of those that matter -- water,
    O2, phosphate, CO2, PPi, ammonia -- appear in over half of all scored
    reactions and several are training anchors. They have to come from
    somewhere, and Zenodo is the only place.
    """
    ours = (row["our_cache_id"] or "").strip()
    has_structure = bool(row["our_smiles"]) and bool(row["our_inchi_key"])
    try:
        n_pkas = int(row["n_pkas"] or 0)
    except ValueError:
        n_pkas = 0
    if not has_structure:
        return "no_structure", "no ModelSEED structure"
    if must_carry is not None:
        if row["seed_id"] in must_carry:
            if not ours:
                return "blocked", "group decomposer refuses it and no cache row exists"
            return "carry", "group decomposer refuses it - Zenodo row is the only source"
        if n_pkas > 0:
            return "create", "build from our SMILES"
        if predicted and row["seed_id"] in predicted:
            # The predictor ran on this structure and returned no ionizable
            # site. That is an answer about the molecule, not missing data, so
            # it is built with an empty pKa list -- terpenes and other
            # hydrocarbons live here. Requires the empty-pKas fix in
            # equilibrator-assets (3e0e4c0); without it creation is refused.
            return "create", "no ionizable groups - build with an empty pKa list"
        return "blocked", "no pKa prediction available"
    if ours:
        return "carry", "cache already holds our structure"
    if n_pkas > 0:
        return "create", "absent from cache - build from our SMILES"
    return "blocked", "absent from cache and no pKas to build it"


def redox_accessions():
    """Accessions referenced by the redox half-reaction training table."""
    import pandas as pd
    from equilibrator_cache.zenodo import get_cached_filepath
    from component_contribution import TRAINING_DATA_REDOX_SETTINGS

    df = pd.read_csv(get_cached_filepath(TRAINING_DATA_REDOX_SETTINGS))
    accs = set()
    for col in ("CID_ox", "CID_red"):
        accs.update(a for a in df[col].dropna() if ":" in str(a))
    return accs


def copy_rows(src, dst, table, ids, id_column="id", chunk=900):
    """Copy whole rows verbatim, preserving primary keys."""
    if not ids:
        return 0
    cols = [r[1] for r in src.execute(f"pragma table_info({table})")]
    collist = ",".join(cols)
    placeholders = ",".join("?" * len(cols))
    ids = list(ids)
    n = 0
    for i in range(0, len(ids), chunk):
        part = ids[i : i + chunk]
        q = (f"select {collist} from {table} where {id_column} in "
             f"({','.join('?' * len(part))})")
        rows = src.execute(q, part).fetchall()
        dst.executemany(
            f"insert or ignore into {table} ({collist}) values ({placeholders})",
            rows)
        n += len(rows)
    dst.commit()
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--max-create", type=int, default=None,
                    help="cap creations, for a smoke run")
    ap.add_argument("--batch", type=int, default=200)
    ap.add_argument("--pkas", type=Path, default=None,
                    help="MolGpKa-format table whose values replace the mapping's "
                         "pKas. Lets caches be built that differ in nothing but "
                         "the pKa source, which is the only way to isolate the "
                         "effect of that source on computed energies.")
    ap.add_argument("--pkas-mode", choices=("nonmarvin", "all"), default="all",
                    help="'all' (default) replaces every row the table covers, "
                         "including rows currently sourced from ChemAxon Marvin. "
                         "That is the arm that measures what an open-source pKa "
                         "layer costs. 'nonmarvin' leaves Marvin rows alone, "
                         "isolating one MolGpKa checkpoint against another -- on "
                         "this mapping that is only 266 compounds, so it answers "
                         "a narrower question than it appears to.")
    ap.add_argument("--no-carry-over", action="store_true",
                    help="build every compound from ModelSEED structures, "
                         "falling back to a Zenodo row only for the allowlist "
                         "in --must-carry. Without this, 22,326 compounds "
                         "import Marvin-derived microspecies from Zenodo.")
    ap.add_argument("--must-carry", type=Path,
                    default=ROOT / "data" / "must_carry_over.tsv",
                    help="compounds the group decomposer refuses, which "
                         "therefore cannot be built; produced by "
                         "tools/must_carry_over.py")
    ap.add_argument("--no-reprotonate", action="store_true",
                    help="leave carried-over rows with Zenodo's ChemAxon-derived "
                         "dissociation constants instead of substituting ours")
    ap.add_argument("--max-create-failures", type=int, default=2000,
                    help="exit non-zero if more CREATIONS fail than this. The "
                         "known-blocked no-pKa population is reported but not "
                         "counted -- it is an input gap, not a build failure.")
    args = ap.parse_args()

    import modelseed_pkas as msp

    plan = list(csv.DictReader(MAPPING.open(), delimiter="\t"))

    predicted = None
    if args.pkas:
        # every compound the predictor was RUN on; those with a structure but
        # no row in its output genuinely have no ionizable group
        predicted = {r["seed_id"] for r in plan
                     if r["our_smiles"] and r["our_inchi_key"]}

    must_carry = None
    if args.no_carry_over:
        if not args.must_carry.exists():
            raise SystemExit(
                f"--no-carry-over needs {args.must_carry}; generate it with "
                f"tools/must_carry_over.py")
        with args.must_carry.open() as fh:
            must_carry = {r["seed_id"] for r in csv.DictReader(fh, delimiter="\t")}
        print(f"no-carry-over mode: {len(must_carry):,} compounds allowlisted "
              f"as un-buildable, everything else built from our structures")

    if args.pkas:
        # Marvin-sourced rows are left alone: substituting there would confound
        # "which MolGpKa model" with "MolGpKa instead of Marvin".
        table = {}
        for row in csv.DictReader(args.pkas.open(), delimiter="\t"):
            vals = []
            for tok in row["value"].split(";"):
                if tok:
                    try:
                        vals.append(float(tok.rsplit(":", 1)[-1]))
                    except ValueError:
                        pass
            table.setdefault(row["external_id"], []).extend(vals)
        swapped = 0
        for row in plan:
            if (args.pkas_mode == "nonmarvin"
                    and (row.get("pka_source") or "").startswith("marvin")):
                continue
            vals = table.get(row["seed_id"])
            if not vals:
                continue
            # -2..16, not 0..14. eQuilibrator counts sites above the
            # reported pH to place the major microspecies, so a group at
            # pKa 14.9 is protonated at pH 7 and carries a proton; dropping
            # it decrements major_ms_index and moves dG'0 by a full
            # RT ln(10) x pH = 9.55 kcal/mol. The old window discarded 210
            # such sites, moving 246 scored reactions, and left 181 compounds
            # with no sites at all -- which then read as "does not ionize",
            # turning a filter choice into a false claim about the molecule.
            # -2..16 is the physical-plausibility bound: nothing in aqueous
            # solution ionizes meaningfully outside it.
            keep = sorted((v for v in vals if -2.0 <= v <= 16.0), reverse=True)
            row["pkas"] = ";".join(f"{v:.4f}" for v in keep)
            row["n_pkas"] = str(len(keep))
            row["pka_source"] = f"molgpka:{args.pkas.stem}"
            swapped += 1
        print(f"pKa override from {args.pkas.name} (mode={args.pkas_mode}): "
              f"replaced {swapped:,} rows")
    buckets = defaultdict(list)
    for row in plan:
        action, reason = classify(row, must_carry, predicted)
        buckets[action].append((row, reason))

    print("plan")
    for k in ("carry", "create", "blocked", "no_structure"):
        print(f"  {k:<14}{len(buckets[k]):7,}")

    # Accessions the training data references; those with no ModelSEED route
    # must be carried over or the training set shrinks silently.
    # build_stage2_cache.training_accessions() reads TECRDB + formation only.
    # The 13 redox couples reference their own metanetx accessions in CID_ox /
    # CID_red, and omitting them makes Trainer.read_redox raise ParseException
    # at train time -- loud, but only after the whole cache has been built.
    # (The "831 training accessions" figure in CACHE_STRATEGY.md inherits the
    # same omission.)
    from tools.build_stage2_cache import training_accessions
    tr_acc = set(training_accessions()) | redox_accessions()
    print(f"  training accessions referenced: {len(tr_acc):,}"
          f"  (incl. redox)")

    if args.dry_run:
        print("\n(dry run - nothing written)")
        return

    if args.out.resolve() == SOURCE_CACHE.resolve():
        sys.exit("refusing to write to the pristine Zenodo cache")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        args.out.unlink()

    from sqlalchemy import create_engine
    from equilibrator_cache.models import Base
    Base.metadata.create_all(create_engine(f"sqlite:///{args.out}"))
    print(f"\ncreated empty cache at {args.out}")

    src = sqlite3.connect(f"file:{SOURCE_CACHE}?mode=ro", uri=True)
    dst = sqlite3.connect(args.out)

    # Registries first: get_compounds looks up the `synonyms` registry and
    # raises NoResultFound against an empty table.
    n = copy_rows(src, dst, "registries", [r[0] for r in src.execute("select id from registries")])
    print(f"  registries copied: {n}")

    # --- carry-over -------------------------------------------------------
    carry_ids = {int(r["our_cache_id"]) for r, _ in buckets["carry"] if r["our_cache_id"]}

    # training-only compounds: resolve each accession in the source cache
    # Namespace lookup is case-folded: the training data carries both "coco:"
    # and "COCO:", and an exact match silently drops the latter.
    reg_by_ns = {ns.lower(): rid
                 for rid, ns in src.execute("select id, namespace from registries")}
    training_ids = set()
    unresolved_training = []
    for acc in sorted(tr_acc):
        ns, _, code = acc.partition(":")
        rid = reg_by_ns.get(ns.lower())
        if rid is None:
            unresolved_training.append((acc, "unknown registry"))
            continue
        hits = [c for (c,) in src.execute(
            "select compound_id from compound_identifiers "
            "where registry_id = ? and accession = ?", (rid, code))]
        if not hits:
            unresolved_training.append((acc, "not in source cache"))
            continue
        training_ids.add(msp.pick_compound(hits, training_ids | carry_ids)
                         if len(hits) > 1 else hits[0])
    training_only = training_ids - carry_ids
    print(f"  carry-over compounds        : {len(carry_ids):,}")
    print(f"  training-only compounds     : {len(training_only):,}"
          f"   (unresolved: {len(unresolved_training)})")

    keep = carry_ids | training_only
    t0 = time.time()
    print(f"  copying {len(keep):,} compounds ...", flush=True)
    copy_rows(src, dst, "compounds", keep)
    for table, col in (("compound_microspecies", "compound_id"),
                       ("compound_identifiers", "compound_id"),
                       ("magnesium_dissociation_constant", "compound_id")):
        try:
            copy_rows(src, dst, table, keep, id_column=col)
        except sqlite3.OperationalError:
            pass
    print(f"    done in {time.time() - t0:.0f}s")

    # Every identifier came from the source cache, including its seed:
    # namespace. Purge it -- accessions are re-derived below from our mapping,
    # never inherited.
    seed_reg = reg_by_ns[SEED_NAMESPACE]
    gone = dst.execute("delete from compound_identifiers where registry_id = ?",
                       (seed_reg,)).rowcount
    dst.commit()
    print(f"  purged inherited seed: identifiers: {gone:,}")
    dst.close()
    src.close()

    # --- creation ---------------------------------------------------------
    to_create = [(r, why) for r, why in buckets["create"]]
    if args.max_create is not None:
        to_create = to_create[: args.max_create]
    by_smiles = defaultdict(list)
    for r, _ in to_create:
        by_smiles[r["our_smiles"]].append(r)
    unique = sorted(by_smiles, key=lambda s: by_smiles[s][0]["seed_id"])
    print(f"\ncreating {len(unique):,} structures for {len(to_create):,} accessions")

    created = {}
    if unique:
        from equilibrator_assets.local_compound_cache import LocalCompoundCache
        batch_casualties = []
        lcc = LocalCompoundCache(pkas_provided_externally=True)
        lcc.load_cache(str(args.out))
        t0 = time.time()
        for i in range(0, len(unique), args.batch):
            chunk = unique[i : i + args.batch]
            pkas = {s: [float(p) for p in by_smiles[s][0]["pkas"].split(";") if p]
                    for s in chunk}
            try:
                results = lcc.get_compounds(chunk, specified_pkas=pkas)
            except Exception as exc:
                # One unparseable structure used to take its whole batch with
                # it -- 199 good compounds lost to one bad cobalamin SMILES.
                # Worse, batch boundaries shift when the input list changes, so
                # two otherwise-identical builds lost *different* compounds and
                # were not comparable. Retry singly and lose only the culprit.
                print(f"  batch {i // args.batch}: FAILED "
                      f"({type(exc).__name__}: {exc}) -- retrying singly",
                      flush=True)
                results = []
                for smiles in chunk:
                    try:
                        results.extend(lcc.get_compounds(
                            [smiles], specified_pkas={smiles: pkas[smiles]}))
                    except Exception as one:
                        batch_casualties.append((smiles, f"{type(one).__name__}: {one}"))
                        results.append(None)
                print(f"    recovered {sum(r is not None for r in results)}"
                      f"/{len(chunk)}", flush=True)
            for smiles, res in zip(chunk, results):
                if res is not None and res.compound is not None:
                    created[smiles] = res.compound.id
            done = min(i + args.batch, len(unique))
            rate = (time.time() - t0) / done
            print(f"  {done}/{len(unique)} structures, {len(created)} created, "
                  f"{rate:.2f}s each, ETA {rate * (len(unique) - done) / 60:.0f} min",
                  flush=True)
        lcc.ccache.session.close()
        if batch_casualties:
            print(f"  {len(batch_casualties)} structure(s) unconvertible:", flush=True)
            for smiles, why in batch_casualties[:10]:
                print(f"    {smiles[:70]}... {why[:90]}", flush=True)

    # --- re-protonate carried-over rows ----------------------------------
    # A carried-over row is Zenodo's, and Zenodo is built by shelling out to
    # ChemAxon cxcalc, so its dissociation constants are Marvin's. That is not
    # forced on us. "The group decomposer refuses this compound" means it has
    # an EMPTY group vector and is handled as a pure reactant-contribution
    # anchor -- not that it cannot be represented. The row supplies id,
    # atom_bag, inchi_key and group_vector=[], all structural facts; only the
    # protonation is ChemAxon's, and we can replace it in place.
    #
    # The microspecies must be regenerated rather than patched: major_ms_index
    # and every ddg_over_rt derive from the pKa list, so editing the constants
    # without rebuilding the ladder leaves the two inconsistent.
    if not args.no_reprotonate:
        dst = sqlite3.connect(args.out)
        moved = kept = emptied = mg_skipped = 0
        for row, _ in buckets["carry"]:
            cid = row.get("our_cache_id")
            if not cid:
                continue
            cid = int(cid)
            ours = [float(x) for x in (row.get("pkas") or "").split(";") if x]
            existing = dst.execute(
                "select dissociation_constants from compounds where id=?", (cid,)
            ).fetchone()
            if existing is None:
                continue
            try:
                theirs = pickle.loads(existing[0]) if isinstance(
                    existing[0], (bytes, bytearray)) else existing[0]
            except Exception:
                theirs = None
            if not theirs:
                emptied += 1          # nothing ChemAxon-derived to displace
                continue
            if not ours:
                kept += 1             # no source of ours; theirs stands, and
                continue              # the provenance audit will report it

            # Magnesium species are indexed by proton count and were computed
            # against the STORED ladder. Our ladders are often shorter --
            # Alberty tabulates only the dissociations that matter between
            # about pH 5 and 9, so phosphate arrives as [7.22] rather than
            # [12.35, 7.20, 2.15] -- and a shorter ladder spans fewer proton
            # counts. Substituting it orphans every Mg species outside the new
            # range, and populate_microspecies then fails outright looking for
            # a reference that no longer exists.
            #
            # So re-protonate only when the new ladder still covers every
            # proton count the Mg data references. Skipping is the honest
            # outcome: a truncated ladder really is less information than the
            # one it would replace, however much better its individual values.
            mg_protons = [
                r[0] for r in dst.execute(
                    "select number_protons from magnesium_dissociation_constant "
                    "where compound_id=?", (cid,))
            ]
            if mg_protons:
                bag = dst.execute(
                    "select atom_bag from compounds where id=?", (cid,)).fetchone()[0]
                try:
                    n_h = (pickle.loads(bag) if isinstance(bag, (bytes, bytearray))
                           else bag).get("H", 0)
                except Exception:
                    n_h = 0
                major = sum(1 for v in ours if v > 7.0)
                reachable = {i - major + n_h for i in range(len(ours) + 1)}
                if not set(mg_protons) <= reachable:
                    mg_skipped += 1
                    continue
            dst.execute(
                "update compounds set dissociation_constants=? where id=?",
                (pickle.dumps(sorted(ours, reverse=True)), cid),
            )
            dst.execute("delete from compound_microspecies where compound_id=?", (cid,))
            moved += 1
        dst.commit()
        dst.close()
        print(f"\nre-protonated carried-over rows")
        print(f"  replaced with our pKas   {moved:,}")
        print(f"  already had none stored  {emptied:,}")
        print(f"  no source of ours, kept  {kept:,}")
        print(f"  ladder too short for Mg  {mg_skipped:,}")
        print(f"  residual ChemAxon exposure: {kept + mg_skipped:,}")
        if moved:
            from equilibrator_assets.thermodynamics import populate_microspecies
            from sqlalchemy import create_engine
            from sqlalchemy.orm import sessionmaker

            eng = create_engine(f"sqlite:///{args.out}")
            with sessionmaker(bind=eng)() as session:
                populate_microspecies(session)
                session.commit()
            eng.dispose()
            print("  microspecies regenerated for the rows that changed")

    # --- attach seed: accessions -----------------------------------------
    dst = sqlite3.connect(args.out)
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    inserts, failures = [], []
    stats = Counter()
    for action in ("carry", "create", "blocked", "no_structure"):
        for r, reason in buckets[action]:
            cid = None
            if action == "carry":
                cid = int(r["our_cache_id"])
            elif action == "create":
                cid = created.get(r["our_smiles"])
                if cid is None:
                    failures.append((r["seed_id"], r["our_inchi_key"],
                                     r["formula"], "creation failed", r["our_smiles"]))
            elif action == "blocked":
                failures.append((r["seed_id"], r["our_inchi_key"], r["formula"],
                                 "no pKas to build it", r["our_smiles"]))
            if cid:
                inserts.append((now, now, cid, seed_reg, r["seed_id"]))
                stats[action] += 1
            else:
                stats[f"{action}: unattached"] += 1
    dst.executemany(
        "insert into compound_identifiers "
        "(created_on, updated_on, compound_id, registry_id, accession) "
        "values (?,?,?,?,?)", inserts)
    dst.commit()
    print(f"\nattached {len(inserts):,} seed: accessions")
    for k, v in sorted(stats.items()):
        print(f"    {k:<28}{v:7,}")

    FAILURES.parent.mkdir(parents=True, exist_ok=True)
    with FAILURES.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["seed_id", "inchi_key", "formula", "reason", "smiles"])
        w.writerows(sorted(failures))
    reasons = Counter(f[3] for f in failures)
    print(f"  wrote {FAILURES}  ({len(failures):,} accessions with no compound)")
    for k, v in reasons.most_common():
        print(f"    {k:<28}{v:7,}")

    # --- repoint kegg: accessions onto our structures ---------------------
    # Attaching seed: accessions is only half the job. The training reactions
    # resolve through kegg:/metanetx: accessions, so where our structure
    # differs the measured reactant contribution stays attached to
    # eQuilibrator's compound and our seed: accession points at one with no
    # measured energy. That is the Path A / Path B difference, and skipping it
    # costs 11 anchors (Nc 631 rather than 642).
    #
    # Only kegg: can be moved -- metanetx.chemical: and coco: have no ModelSEED
    # alias route.
    import re as _re
    structures = msp.load_structures()
    kegg2seed = {}
    for cpd, rec in structures.items():
        for al in rec["aliases"]:
            if _re.fullmatch(r"C\d{5}", al):
                kegg2seed.setdefault(al, []).append(cpd)
    kegg_reg = reg_by_ns["kegg"]
    seed_target = {a: c for a, c in dst.execute(
        "select accession, compound_id from compound_identifiers where registry_id=?",
        (seed_reg,))}
    kegg_target = {a: c for a, c in dst.execute(
        "select accession, compound_id from compound_identifiers where registry_id=?",
        (kegg_reg,))}
    moves = []
    rp = Counter()
    for acc in sorted(tr_acc):
        if not acc.lower().startswith("kegg:"):
            rp["not a kegg accession"] += 1
            continue
        kid = acc.split(":", 1)[1]
        cur = kegg_target.get(kid)
        if cur is None:
            rp["kegg accession absent"] += 1
            continue
        targets = {seed_target[s] for s in kegg2seed.get(kid, []) if s in seed_target}
        if not targets:
            rp["no ModelSEED compound"] += 1
            continue
        if len(targets) > 1:
            rp["ambiguous"] += 1
            continue
        tgt = targets.pop()
        if tgt == cur:
            rp["already agrees"] += 1
            continue
        rp["REPOINTED to our structure"] += 1
        moves.append((tgt, kegg_reg, kid))
    dst.executemany("update compound_identifiers set compound_id=? "
                    "where registry_id=? and accession=?", moves)
    dst.commit()
    print(f"\nrepointing kegg: accessions onto our structures")
    for k, v in rp.most_common():
        print(f"    {k:<32}{v:6,}")

    # --- verification ------------------------------------------------------
    print("\nverifying")
    dup = dst.execute(
        "select count(*) from (select accession from compound_identifiers "
        "where registry_id = ? group by accession having count(*) > 1)",
        (seed_reg,)).fetchone()[0]
    print(f"  seed: accessions mapping to >1 compound : {dup}")
    ncpd = dst.execute("select count(*) from compounds").fetchone()[0]
    nseed = dst.execute("select count(*) from compound_identifiers where registry_id = ?",
                        (seed_reg,)).fetchone()[0]
    print(f"  compounds in cache                      : {ncpd:,}")
    print(f"  seed: identifiers                       : {nseed:,}")
    for name, cpd in (("water", "cpd00001"), ("phosphate", "cpd00009"),
                      ("CO2", "cpd00011"), ("ammonia", "cpd00013"),
                      ("pyruvate", "cpd00020")):
        hit = dst.execute(
            "select compound_id from compound_identifiers "
            "where registry_id = ? and accession = ?", (seed_reg, cpd)).fetchone()
        print(f"    {name:<10}{cpd}  -> {hit[0] if hit else 'MISSING'}")
    dst.close()

    if dup:
        sys.exit(f"FATAL: {dup} seed: accessions map to more than one compound")
    # Only unexpected failures gate the run. The no-pKa population is a known
    # input gap (6,405 ModelSEED compounds have never been through MolGpKa);
    # counting it here would make every run fail for a reason the run cannot
    # fix, which trains people to ignore the exit code.
    n_create_fail = reasons.get("creation failed", 0)
    if n_create_fail > args.max_create_failures:
        sys.exit(f"FATAL: {n_create_fail} structures failed to build, over the "
                 f"--max-create-failures limit of {args.max_create_failures}. "
                 f"See {FAILURES}.")


if __name__ == "__main__":
    main()
