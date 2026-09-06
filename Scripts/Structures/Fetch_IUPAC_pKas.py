#!/usr/bin/env python
"""Build the local IUPAC pKa flat file. Its output is deliberately gitignored.

The IUPAC Digitized pKa Dataset (Zenodo 10.5281/zenodo.15375522) is 24,211
measured pKa values digitized with IUPAC's permission from Serjeant & Dempsey,
*Ionisation Constants of Organic Acids in Aqueous Solution* (1979) and Perrin,
*Dissociation Constants of Organic Bases in Aqueous Solution* (1965, with the
1972 supplement). Every row carries SMILES and InChI, so compounds are matched
structurally rather than by name, and values are typed ``pKa1`` / ``pKa2`` /
``pKaH1``, which makes the macroscopic ladder explicit at source -- exactly what
eQuilibrator's transform requires and what a per-site predictor cannot supply.

**The dataset is CC-BY-NC-4.0 and is never committed.** This database is
redistributed without a non-commercial restriction, so the values are consumed
by the pipeline and left on local disk. The output path is gitignored; this
script and the ``sources.yaml`` entry are the reproducible record.

Selection, where a compound has several measurements of the same dissociation:

* rows with a cosolvent are dropped -- those are not aqueous values
* 25 C is preferred, then 20 C, then anything with a stated temperature
* among the survivors the **median** is taken, and the spread is reported so a
  consumer can see where the literature disagrees

Ionic strength is the known weakness. eQuilibrator wants pKas at **zero** ionic
strength, applying its own Debye-Huckel correction from the microspecies charge
(``_legendre_transform`` in ``equilibrator_cache.thermodynamic_constants``), so
supplying pre-corrected values would double-count. This dataset has **no ionic
strength column** -- conditions appear only in free-text ``remarks`` such as
``"C=0.04, f+/- taken equal f+/-(KCl)"``. Rows are therefore emitted with the
remark text preserved and an ``ionic_strength_stated`` flag, and the caller
decides. Alberty's BasicBiochemData is uniformly I=0 and needs no such care,
which is why it remains the preferred source where the two overlap.
"""
import argparse
import csv
import re
import statistics as st
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Biochemistry" / "Structures" / "ModelSEED" / "pkas" / "iupac_v2_3b.tsv"
ZENODO = ("https://zenodo.org/records/15375522/files/"
          "IUPAC/Dissociation-Constants-v2-3b.zip?download=1")
MEMBER = "iupac_high-confidence_v2_3.csv"
# free-text markers that a row states an ionic strength or supporting electrolyte
IONIC_HINT = re.compile(r"\bI\s*=|\bKCl\b|\bNaCl\b|\bKNO3\b|ionic", re.I)


def fetch(cache: Path):
    if cache.exists():
        return cache
    cache.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {ZENODO}", flush=True)
    with urlopen(ZENODO) as r, cache.open("wb") as fh:
        fh.write(r.read())
    return cache


def rows_from(zip_path: Path):
    with zipfile.ZipFile(zip_path) as z:
        name = next(n for n in z.namelist() if n.endswith(MEMBER))
        with z.open(name) as fh:
            text = fh.read().decode("utf-8", errors="ignore")
    return list(csv.DictReader(text.splitlines()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", type=Path,
                    default=Path.home() / ".cache" / "iupac_pka.zip")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()

    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    rows = rows_from(fetch(a.zip))
    print(f"IUPAC rows: {len(rows):,}")

    # ModelSEED structures, keyed on the InChIKey connectivity block so that a
    # difference in stored protonation does not prevent a match
    seed = defaultdict(list)
    struct = ROOT / "Biochemistry" / "Structures" / "Unique_ModelSEED_Structures.txt"
    with struct.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Type"] == "InChIKey":
                seed[r["Structure"].split("-")[0]].append(r["ID"])
    print(f"ModelSEED InChIKey skeletons: {len(seed):,}")

    grouped = defaultdict(list)
    dropped = {"cosolvent": 0, "unparseable": 0, "no match": 0, "no value": 0}
    for r in rows:
        if r["cosolvent"].strip():
            dropped["cosolvent"] += 1
            continue
        try:
            value = float(r["pka_value"])
        except (TypeError, ValueError):
            dropped["no value"] += 1
            continue
        mol = Chem.MolFromInchi(r["InChI"]) if r["InChI"] else None
        if mol is None:
            dropped["unparseable"] += 1
            continue
        skel = Chem.MolToInchiKey(mol).split("-")[0]
        hits = seed.get(skel)
        if not hits:
            dropped["no match"] += 1
            continue
        for cpd in hits:
            grouped[(cpd, r["pka_type"])].append((value, r))

    def pick(entries):
        """25 C if any row states it, else 20 C, else whatever has a temperature."""
        for want in ("25", "20"):
            sub = [e for e in entries if e[1]["T"] == want]
            if sub:
                return sub, want
        sub = [e for e in entries if e[1]["T"] not in ("", "not_stated")]
        return (sub or entries), "any"

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["seed_id", "pka_type", "pka_value", "n_measurements",
                    "spread", "temperature_C", "ionic_strength_stated",
                    "assessment", "refs", "remarks"])
        n = 0
        for (cpd, kind), entries in sorted(grouped.items()):
            sub, temp = pick(entries)
            vals = [v for v, _ in sub]
            spread = (max(vals) - min(vals)) if len(vals) > 1 else 0.0
            ionic = any(IONIC_HINT.search(e[1]["remarks"] or "") for e in sub)
            w.writerow([
                cpd, kind, f"{st.median(vals):.3f}", len(vals), f"{spread:.3f}",
                temp, int(ionic),
                ";".join(sorted({e[1]["assessment"] for e in sub if e[1]["assessment"]})),
                ";".join(sorted({e[1]["ref"] for e in sub if e[1]["ref"]}))[:80],
                (sub[0][1]["remarks"] or "")[:100],
            ])
            n += 1

    print(f"\nwrote {a.out}  ({n:,} compound-dissociation rows)")
    print(f"  distinct ModelSEED compounds: {len({c for c, _ in grouped}):,}")
    for k, v in dropped.items():
        print(f"  dropped, {k:<14}{v:>8,}")
    print("\nThis file is gitignored. It is licensed CC-BY-NC-4.0 and must not "
          "be committed or redistributed.")


if __name__ == "__main__":
    main()
