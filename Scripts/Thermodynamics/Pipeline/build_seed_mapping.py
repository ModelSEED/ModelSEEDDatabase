"""Build the seed:cpd##### -> structure -> cache-compound mapping table.

This is the Path A artifact: for every ModelSEED compound, what structure we
hold for it, whether the eQuilibrator cache already contains that structure,
and how the cache's *current* (MetaNetX-derived) seed: mapping compares.

The classification column drives the rebuild:

    exact          our structure is already in the cache under this key, and
                   the existing seed: mapping agrees -- nothing to do
    remap          the cache has our structure, but the existing seed: mapping
                   points somewhere else -- repoint the accession
    stereo         existing mapping differs only in stereo/protonation layer
                   -- needs a per-compound curation call
    wrong_molecule existing mapping is a different molecule entirely
    absent         our structure is not in the cache at all -- create it
    unmapped       no seed: accession in the cache today (pure addition)
    no_pkas        we have a structure but no usable pKas
    no_structure   ModelSEED has no structure for this compound

Writes a TSV and prints a summary. Read-only with respect to the cache.
"""

import csv
import os
import sqlite3
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import modelseed_pkas as msp  # noqa: E402

# The eQuilibrator working tree: caches and fitted parameters, gigabytes, not
# in this repository. Named by environment variable since relocation.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))


CACHE = Path(os.path.expanduser("~/.cache/equilibrator/compounds.sqlite"))
OUT = ROOT / "data" / "seed_mapping.tsv"


def connectivity(key):
    """First InChIKey block: the molecular skeleton."""
    return key.split("-")[0] if key else None


def skeleton_and_stereo(key):
    """First two InChIKey blocks: skeleton plus stereo/isotope layer."""
    return "-".join(key.split("-")[:2]) if key else None


_MOBILE_H = re.compile(r"\(H(\d*)([-+]?),([\d,]+)\)")


def inchi_layers(inchi):
    """formula, /c connectivity, /h hydrogens, /t stereo -- or None."""
    if not inchi or not inchi.startswith("InChI="):
        return None
    parts = inchi.split("/")
    pick = lambda p: next((x for x in parts if x.startswith(p)), "")
    return (parts[1] if len(parts) > 1 else "",
            pick("c"), pick("h"), pick("t"))


def hydrogen_content(h_layer):
    """(fixed spec, mobile proton count, atoms those protons may occupy).

    Two InChIs can describe the same molecule and still hash differently in the
    first InChIKey block, because InChI may partition the mobile hydrogens into
    different groups. FAD is the case: our layer separates the exchangeable
    protons per functional group where the cache merges all eleven into one.
    Comparing content rather than text tolerates that regrouping while still
    refusing a genuine tautomer, where the protons sit somewhere else.
    """
    if not h_layer:
        return "", 0, frozenset()
    i = h_layer.find("(")
    fixed = (h_layer[:i] if i >= 0 else h_layer).rstrip(",")
    total, atoms = 0, set()
    for count, _sign, alist in _MOBILE_H.findall(h_layer[i:] if i >= 0 else ""):
        total += int(count) if count else 1
        atoms.update(int(a) for a in alist.split(",") if a)
    return fixed, total, frozenset(atoms)


def cache_structure_index(sqlite_path):
    """(formula, /c) -> [(compound_id, layers)] for every cache row with an InChI."""
    index = {}
    con = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
    for cid, inchi in con.execute(
            "select id, inchi from compounds where inchi is not null"):
        layers = inchi_layers(inchi)
        if layers:
            index.setdefault((layers[0], layers[1]), []).append((cid, layers))
    con.close()
    return index


def regrouped_match(our_inchi, index):
    """Find a cache row that is the same molecule under a different /h grouping.

    Requires identical formula, connectivity and stereo, an identical fixed-H
    spec, and the same mobile protons over the same atoms. Candidates are
    SEARCHED, never taken first-come: IPP and DMAPP are constitutional isomers
    that InChI separates only in /h, the cache holds a row for each, and
    first-hit selection would hand DMAPP the IPP row.
    """
    layers = inchi_layers(our_inchi)
    if not layers:
        return None
    ours = hydrogen_content(layers[2])
    for cid, cand in index.get((layers[0], layers[1]), []):
        if cand[3] != layers[3]:
            continue                      # stereo must agree
        if hydrogen_content(cand[2]) == ours:
            return cid
    return None


def compare_keys(a, b):
    """Classify how two InChIKeys differ.

    An InChIKey is SKELETON-STEREO-P. The final single character encodes
    protonation state, so a difference confined to it means the same chemical
    entity at a different proton count -- which is expected here, since
    eQuilibrator stores whichever microspecies its pipeline produced while
    ModelSEED stores the Marvin pH-7 form. Those are not mapping errors.
    """
    if a == b:
        return "same"
    if connectivity(a) != connectivity(b):
        return "different_skeleton"
    if skeleton_and_stereo(a) != skeleton_and_stereo(b):
        return "different_stereo"
    return "different_protonation"


def main():
    print("loading ModelSEED structures + pKas from pinned dev export ...")
    structures = msp.load_structures()
    # Deliberately WITHOUT the "cache" tier. That tier resolves ladders through
    # our_cache_id, which is what this script writes, so including it here would
    # make the mapping depend on its own previous output -- two consecutive runs
    # could disagree. The mapping is defined by structures plus our own pKa
    # sources; the cache tier is layered on afterwards by write_resolved_pkas.
    dataset = msp.build_dataset(
        preference=tuple(s for s in msp.DEFAULT_PKA_PREFERENCE if s != "cache"))

    print("indexing cache structures for regrouped-mobile-H matching ...")
    struct_index = cache_structure_index(CACHE)

    print("reading compound cache ...")
    training = msp.training_compound_ids()
    by_key = msp.cache_inchi_keys(CACHE, training)
    # Same preference rule for the connectivity-only fallback: a training
    # compound wins, otherwise the lowest id. Plain setdefault would hand back
    # whichever happened to be encountered first.
    # Match on skeleton+stereo, i.e. ignoring only the protonation block.
    # Connectivity-only matching would accept a different stereoisomer as "our
    # structure", which under a ModelSEED-authoritative policy is exactly the
    # substitution we are trying to stop making.
    ss_groups = {}
    for ik, cid in by_key.items():
        ss_groups.setdefault(skeleton_and_stereo(ik), []).append(cid)
    by_conn = {c: msp.pick_compound(v, training) for c, v in ss_groups.items()}
    current_seed = msp.cache_seed_identifiers(CACHE)
    cache_formula = msp.cache_formulas(CACHE)

    rows = []
    counts = {}
    for cpd in sorted(structures):
        rec = structures[cpd]
        our_key = rec["inchi_key"]
        our_smiles = rec["smiles"]
        entry = dataset.get(cpd)

        cur = current_seed.get(cpd)
        cur_id, cur_key = (cur if cur else (None, None))

        # Where does OUR structure live in the cache?
        our_id = by_key.get(our_key) if our_key else None
        if our_id is None and our_key:
            our_id = by_conn.get(skeleton_and_stereo(our_key))
            matched_on = "protonation-insensitive" if our_id else ""
        else:
            matched_on = "exact" if our_id else ""
        # Last resort: the same molecule whose InChI groups its mobile
        # hydrogens differently, which changes the first InChIKey block and so
        # defeats both matches above. FAD and FMN are the notable cases.
        if our_id is None and rec.get("inchi"):
            our_id = regrouped_match(rec["inchi"], struct_index)
            if our_id is not None:
                matched_on = "mobile-h-regrouped"

        diff = compare_keys(cur_key, our_key) if cur_key else None

        if not our_smiles or not our_key:
            cls = "no_structure"
        elif entry is None:
            cls = "no_pkas"
        elif cur is None:
            cls = "unmapped_present" if our_id is not None else "absent"
        elif diff == "same":
            cls = "exact"
        elif diff == "different_protonation":
            # Same compound, different proton count. eQuilibrator models
            # protonation via microspecies, so the accession still points at
            # the right chemical entity -- benign.
            cls = "protonation_only"
        elif diff == "different_stereo":
            cls = "stereo_conflict"
        elif matched_on == "mobile-h-regrouped":
            # Same molecule, different mobile-H partition. Verified on content,
            # not merely on the hash disagreeing.
            cls = "regrouped_mobile_h"
        elif our_id is not None:
            cls = "remap"
        else:
            cls = "wrong_molecule"

        # For skeleton-level disagreements, does the formula agree? Same
        # formula means isomer / tautomer / ring form, not a wrong mapping.
        formula_match = ""
        if cls in ("remap", "wrong_molecule") and cur_id is not None:
            cf = cache_formula.get(int(cur_id))
            of = msp.normalize_formula(rec["formula"])
            if cf and of:
                formula_match = "same" if cf == of else "different"

        counts[cls] = counts.get(cls, 0) + 1
        rows.append(
            {
                "seed_id": cpd,
                "classification": cls,
                "our_inchi_key": our_key or "",
                "our_smiles": our_smiles or "",
                "our_cache_id": our_id if our_id is not None else "",
                "matched_on": matched_on,
                "our_id_is_training": int(our_id in training) if our_id else 0,
                "current_cache_id": cur_id if cur_id is not None else "",
                "current_inchi_key": cur_key or "",
                "current_is_training": int(cur_id in training) if cur_id else 0,
                "key_diff": diff or "",
                "formula_match": formula_match,
                "n_pkas": len(entry["pkas"]) if entry else 0,
                "pka_source": entry["pka_source"] if entry else "",
                "pkas": ";".join(str(p) for p in entry["pkas"]) if entry else "",
                "formula": rec["formula"],
                "charge": rec["charge"],
            }
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nwrote {len(rows)} rows to {OUT}\n")
    width = max(len(k) for k in counts)
    for cls, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {cls:<{width}}  {n:6d}")

    # The subset that actually needs work, and the risk cases.
    real_conflicts = ("remap", "stereo_conflict", "wrong_molecule")
    touching_training = [
        r for r in rows
        if r["current_is_training"] and r["classification"] in real_conflicts
    ]
    benign_training = [
        r for r in rows
        if r["current_is_training"] and r["classification"] == "protonation_only"
    ]
    print(f"\n  training compounds, protonation-only difference (benign): "
          f"{len(benign_training)}")
    print(f"  training compounds needing a human call: {len(touching_training)}")
    for r in touching_training:
        print(f"    {r['seed_id']:<10} cache_id={r['current_cache_id']:<7} "
              f"{r['classification']:<15} {r['current_inchi_key']:<29} "
              f"vs {r['our_inchi_key']}")


if __name__ == "__main__":
    main()
