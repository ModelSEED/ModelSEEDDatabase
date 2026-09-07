#!/usr/bin/env python
"""Emit the per-compound protonation provenance that Figure 2A reports.

WHY THIS EXISTS. The paper states what fraction of the protonation layer is
open-source. Until 2026-09-07 those percentages were CONSTANTS TYPED INTO
make_figures.py, derived from a file that is not shipped
(eQuilibrator/data/resolved_pkas.provenance.tsv). A reader could not reproduce
them from anything in this repository, which is not acceptable in a paper whose
argument is transparency. This writes one small table that makes both the
by-compound and the by-reaction figures reproducible from released files.

THE CLASSIFICATION, and why it is not simply the pKa cascade's answer.

A compound reaches the shipped cache one of two ways:

  BUILT      -- we created the row from a ModelSEED structure, and the pKa
                cascade chose a tier: alberty > iupac > cache > marvin >
                molgpka, with the two ChemAxon-derived tiers GATED so they are
                consulted only where MolGpKa provably degenerates.

  CARRIED OVER -- the pinned eQuilibrator release already held a row with a
                matching structure, so the builder reused it wholesale:
                microspecies, pKas and group vector. The cascade never ran.

That second path is why the old accounting was wrong. It counted the cascade's
answers and treated carry-over as a minor category (525 compounds); the cache
actually holds 4,122 rows dated to the pinned release. Those compounds carry
ChemAxon-derived ladders no matter what our predictor would have said, so they
belong in the ChemAxon column. Counting them correctly moves the layer from a
reported 20% ChemAxon-derived to 28%, and from 87% of scored reactions to 95%.
The dependency is larger than the paper claimed, not smaller.

INPUTS
    resolved_pkas.provenance.tsv   the cascade's per-compound choice
                                   (PKA_PROVENANCE, default in $EQUILIBRATOR_DIR)
    cache_final/compounds.sqlite   for the created/carried-over split
                                   (row date: the pinned release is 2019)

OUTPUT
    Biochemistry/Thermodynamics/ProtonationEvidence/pka_provenance.tsv
Tally its `shipped_provenance` column to reproduce Figure 2A's left bar; join
it to ModelSEED_Reaction_Energies.tsv on the compounds in `formula` to
reproduce the right bar.
"""
import argparse, csv, os, sqlite3, sys
from pathlib import Path

HERE = Path(__file__).resolve()
REPO = HERE.parents[3]
EQ = Path(os.environ.get("EQUILIBRATOR_DIR", "/scratch/seaver/Claude_Projects/eQuilibrator"))
CHEMAXON = {"carried_over", "marvin", "cache"}   # all three descend from cxcalc
OPEN_SRC = {"molgpka", "iupac", "alberty"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", type=Path, default=EQ / "data/cache_final/compounds.sqlite")
    ap.add_argument("--pka-provenance", type=Path,
                    default=Path(os.environ.get("PKA_PROVENANCE",
                                 EQ / "data/resolved_pkas.provenance.tsv")))
    ap.add_argument("--out", type=Path,
                    default=REPO / "Biochemistry/Thermodynamics/ProtonationEvidence/pka_provenance.tsv")
    a = ap.parse_args()
    for p, what in ((a.cache, "cache"), (a.pka_provenance, "pKa provenance table")):
        if not p.exists():
            sys.exit(f"missing {what}: {p}")

    sys.path.insert(0, str(EQ / "tools"))
    from modelseed_pkas import cache_seed_identifiers

    seeds = cache_seed_identifiers(str(a.cache))
    con = sqlite3.connect(str(a.cache))
    origin = {cid: co[:4] for cid, co in con.execute("select id, created_on from compounds")}
    resolved = {}
    with a.pka_provenance.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            resolved[r["seed_id"]] = (r["source"], r.get("n_pkas", ""))

    rows = []
    for seed, (cid, _k) in sorted(seeds.items()):
        carried = origin.get(cid) == "2019"
        src, npk = resolved.get(seed, ("unresolved", ""))
        # A carried-over row keeps the pinned release's ladder whatever the
        # cascade would have chosen, so carry-over wins the classification.
        eff = "carried_over" if carried else src
        rows.append({
            "seed_id": seed, "cache_id": cid,
            "cache_row": "carried_over" if carried else "built",
            "pka_cascade_source": src,
            "effective_source": eff,
            "shipped_provenance": ("chemaxon" if eff in CHEMAXON
                                   else "open" if eff in OPEN_SRC else "unresolved"),
            "n_pkas": npk})

    a.out.parent.mkdir(parents=True, exist_ok=True)
    import collections
    tal = collections.Counter(r["shipped_provenance"] for r in rows)
    eff = collections.Counter(r["effective_source"] for r in rows)
    with a.out.open("w", newline="") as fh:
        fh.write("# Per-compound protonation provenance for the shipped cache.\n"
                 "# shipped_provenance: chemaxon = the ladder descends from ChemAxon cxcalc,\n"
                 "#   whether via a carried-over row from the pinned release or via the gated\n"
                 "#   cache/marvin tiers. open = MolGpKa, IUPAC or Alberty.\n"
                 "# cache_row=carried_over means the builder reused the pinned release's row\n"
                 "#   wholesale and the pKa cascade never ran for that compound.\n"
                 f"# totals: {dict(tal)}\n"
                 f"# by effective source: {dict(eff)}\n")
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    n = len(rows)
    print(f"wrote {a.out}  ({n} compounds)")
    for k, v in tal.most_common():
        print(f"  {k:12s} {v:6d}  {v/n:5.1%}")


if __name__ == "__main__":
    main()
