"""Adapt ModelSEED structures and pKas (as they exist on dev) for eQuilibrator.

Reads the pinned export produced by ``tools/extract_dev_structures.sh``, which
mirrors the layout on ``upstream/dev`` of the ModelSEED biochemistry database:

    Unique_ModelSEED_Structures.txt        cpd-keyed, pH-7 charged structures
    <SOURCE>/pkas/marvin_23.4.tsv          Marvin pKa/pKb, keyed by source id
    <SOURCE>/protonations/marvin_*_ph7.tsv per-source pH-7 structures
    ModelSEED/pkas/opam2_molgpka.tsv       MolGpKa pKa/pKb, keyed by cpd id

Note this is *not* the layout used by the in-flight
``pubchem-stereo-loss-guard_20260704`` branch, which stores the same data as
``<SOURCE>/pKa_Strings.txt`` and ``SMILE_ChargedStrings.txt``. dev is the
authority; the branch files are not merged.

Structures are taken from the cpd-level aggregate rather than the per-source
files because eQuilibrator needs one structure per ModelSEED compound, already
protonated at pH 7 -- the bypass path derives ``atom_bag`` from whatever
structure it is given and treats its hydrogen count as the major microspecies.
"""

from __future__ import annotations

import os
import pickle
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Set, Tuple
import sys
import csv
import math
import re

# The pinned ModelSEED structure/pKa export the pipeline reads. Lives in the
# eQuilibrator working tree, which is named by environment variable because it
# is not part of this repository.
EQ_DIR = Path(os.environ.get("EQUILIBRATOR_DIR",
                             "/scratch/seaver/Claude_Projects/eQuilibrator"))
DEV_DIR = EQ_DIR / "data" / "msd_dev"
# The upstream (Zenodo) compound cache. build_seed_mapping.py resolves
# ``our_cache_id`` against this database, so anything reading that id must read
# it here too.
UPSTREAM_CACHE = Path(os.path.expanduser("~/.cache/equilibrator/compounds.sqlite"))
MARVIN_SOURCES = ("KEGG", "ChEBI", "MetaCyc", "Rhea")

# Marvin is what eQuilibrator itself shells out to (cxcalc *is* Marvin), so
# Marvin pKas are the like-for-like substitute for the values the upstream
# pipeline would have computed. MolGpKa is a different predictor and gives
# materially different numbers, so it is only used where Marvin has nothing.
# Resolution order, set by measurement rather than preference. Scored against
# Alberty's I=0 ladders on the 15 compounds where the sources overlap, the mean
# |error| is: IUPAC 0.46, MolGpKa 1.35 pKa units (IUPAC closer on 11 of 15).
# Alberty itself is the reference, so it leads.
#
#   alberty  macroscopic, uniformly at zero ionic strength -- which is what
#            eQuilibrator wants, since it applies its own Debye-Huckel
#            correction from the microspecies charge. Ingests with no
#            conversion. Narrow: 41 ModelSEED compounds, but they include
#            ATP, Pi, ADP, CoA and PPi.
#   iupac    macroscopic and measured, 1,453 compounds. Conditions vary and the
#            dataset has no ionic-strength column, so values carry a bounded
#            offset -- ATP reads 6.53 at I=0.1 against Alberty's 7.60 at I=0.
#            Still ~3x closer than a per-site predictor.
#   molgpka  microscopic per-site values on the neutral molecule. Fine where
#            the ionizable sites are chemically distinct; collapses to a single
#            repeated value where they are equivalent, so sugar phosphates are
#            wrong by 4 pKa units.
#
# "cache"    the macroscopic ladder already held in the upstream eQuilibrator
#            cache for the same molecule. ChemAxon-derived, so not open, but it
#            is a SEQUENTIAL ladder rather than per-site values, which is what
#            the transform actually needs: it reproduces measured literature
#            ladders about 79% of the time against MolGpKa's 35%, and the gap is
#            entirely on the polyprotic acids that dominate the database.
#
# The Marvin *tier* is still absent -- we never invoke ChemAxon. The cache tier
# reuses a value already published in a pinned, citable release, which is a
# reproducibility claim ("reproducible by citation") rather than an open-source
# one, and the methods section has to say so. Placing it above molgpka rather
# than below is the 2026-09-03 decision: on symmetric polyprotic sites a
# per-site predictor structurally cannot express the second dissociation, so
# preferring it there would be choosing the worse answer on principle.
DEFAULT_PKA_PREFERENCE = ("alberty", "iupac", "cache", "marvin", "molgpka")

# Tiers that are ChemAxon-derived, and so are consulted ONLY where MolGpKa
# structurally cannot answer. MolGpKa's failure is specific and detectable, not
# general: it predicts once per ionizable atom on the neutral molecule, so
# chemically equivalent sites return the same number and every dissociation
# after the first is a copy. Where the sites are distinct it scores MAE 0.55
# against measurement and there is no reason to displace it.
#
# Gating matters beyond provenance. Substituting a source for one compound
# changes its site count, and if its reaction partner keeps the old
# enumeration the two no longer cancel -- that is what moved rxn08789 by
# 370 kcal/mol when the cache tier was applied unconditionally. Every needless
# substitution is a chance to break a reaction, so the tier is confined to the
# 6,664 compounds (21%) that provably need it.
GATED_TIERS = frozenset({"cache", "marvin"})

# Kept so the pre-2026-09 behaviour can still be reproduced for comparison.
LEGACY_PKA_PREFERENCE = ("marvin", "molgpka")


# A single atom, optionally charged: Mn2+, Zn2+, K+, elemental S. Formula is one
# element symbol with no count, or a count of 1.
_MONATOMIC = re.compile(r"^[A-Z][a-z]?1?$")


def is_bare_ion(rec: dict) -> bool:
    """Does this structure have no ionizable proton by construction?

    A monatomic species cannot donate a proton, so its pKa list is EMPTY -- and
    that is an answer, not missing data. ChemAxon nonetheless returns 3.09 for
    the metal cations; equilibrator-assets calls this out in its own
    COMPOUND_EXCEPTIONS ("Metal Cations get multiple pKa values from ChemAxon,
    which is obviously a bug") and overrides a handful of species, but not the
    wider set.

    This matters specifically because the fallback tiers are gated on "MolGpKa
    said nothing", and MolGpKa correctly says nothing about a bare ion. Without
    this check the gate opens and we would emit the artefact as our own
    resolved value rather than merely inheriting it.
    """
    return bool(_MONATOMIC.match((rec.get("formula") or "").strip()))


REPO = Path(__file__).resolve().parents[3]
EXCLUSIONS = (REPO / "Biochemistry" / "Thermodynamics" /
              "ProtonationEvidence" / "substitution_exclusions.tsv")


def load_substitution_exclusions(path: Path = None) -> Set[str]:
    """Compounds that keep MolGpKa even where a macroscopic ladder exists.

    Large glycans -- LPS and lipid A, peptidoglycan, lipid-linked
    oligosaccharides -- whose reaction partners are equally large and cannot be
    substituted alongside them. Moving one side alone leaves dozens of sugar
    hydroxyls counted on the other, and each unbalanced proton is 9.55 kcal/mol.

    Fitted, not derived: regenerate with find_substitution_exclusions.py after
    any change to the substitution policy or the pKa sources.
    """
    path = path or EXCLUSIONS
    if not path.exists():
        return set()
    out = set()
    with path.open() as fh:
        for row in csv.DictReader((l for l in fh if not l.startswith("#")),
                                  delimiter="\t"):
            if row.get("seed_id"):
                out.add(row["seed_id"])
    return out


def collapses(p_kas: Iterable[float], places: int = 2) -> bool:
    """Does this look like MolGpKa's symmetry collapse?

    A repeated value means two chemically equivalent sites returned the same
    prediction, which is the signature of the model never seeing the charge
    left behind by the previous deprotonation. Orthophosphate comes back as
    2.11 three times where the true ladder is 2.15 / 7.20 / 12.35.
    """
    rounded = [round(float(v), places) for v in p_kas]
    return len(rounded) != len(set(rounded))


def needs_macroscopic_fallback(cpd: str, molgpka: Dict[str, List[float]]) -> bool:
    """True when MolGpKa cannot answer for this compound.

    Either it collapses, or it has nothing to say at all -- in which case a
    fallback is not displacing anything.
    """
    values = molgpka.get(cpd)
    if not values:
        return True
    return collapses(values)


def normalize_pkas(
    p_kas: Iterable[float], min_ph: float = -2.0, max_ph: float = 16.0
) -> List[float]:
    """Filter to (min_ph, max_ph) and sort descending.

    Mirrors ``equilibrator_assets.generate_compound._normalize_pkas`` so the
    values handed over are already in the order ``create_microspecies_mappings``
    assumes when it slices the list.
    """
    return sorted((p for p in map(float, p_kas) if min_ph < p < max_ph), reverse=True)


def _parse_pka_value(field: str) -> List[float]:
    """Pull the values out of a ``frag:atom:value;frag:atom:value`` field."""
    values = []
    for triple in field.strip().split(";"):
        if not triple:
            continue
        try:
            values.append(float(triple.rsplit(":", 1)[-1]))
        except ValueError:
            continue
    return values


def _read_pka_table(path: Path) -> Dict[str, List[float]]:
    """Read one pKa TSV, merging the pKa and pKb rows for each id.

    cxcalc is invoked upstream with both ``--na`` and ``--nb`` and its acidic
    and basic columns are concatenated into a single dissociation_constants
    list, so merging the two row kinds reproduces that.
    """
    merged: Dict[str, List[float]] = defaultdict(list)
    if not path.is_file():
        return {}
    with path.open() as handle:
        header = next(handle, "").rstrip("\n").split("\t")
        try:
            i_id, i_kind, i_val = (
                header.index("external_id"),
                header.index("kind"),
                header.index("value"),
            )
        except ValueError:
            return {}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_id, i_kind, i_val):
                continue
            if parts[i_kind] not in ("pKa", "pKb"):
                continue
            merged[parts[i_id]].extend(_parse_pka_value(parts[i_val]))
    return dict(merged)


def load_marvin_pkas(dev_dir: Path = DEV_DIR) -> Dict[str, Dict[str, List[float]]]:
    """Marvin pKas per source, keyed by that source's own external id."""
    return {
        source: _read_pka_table(Path(dev_dir) / source / "pkas" / "marvin_23.4.tsv")
        for source in MARVIN_SOURCES
    }


def load_molgpka_pkas(dev_dir: Path = DEV_DIR) -> Dict[str, List[float]]:
    """MolGpKa pKas, already keyed by ModelSEED cpd id.

    The file is whichever ModelSEED pKa table ``sources.yaml`` flags
    ``consumed_by_production``, so a rerun is a manifest edit rather than a
    code edit. Falls back to the original opam2 export if the manifest cannot
    be read -- but note that export covers 22,395 compounds against the
    2026-09 run's 32,201, and records the ionizable HYDROGEN rather than the
    heavy atom bearing it, so silently using it would understate coverage by a
    third.
    """
    base = Path(dev_dir) / "ModelSEED" / "pkas"
    manifest = Path(dev_dir) / "sources.yaml"
    chosen = None
    if manifest.exists():
        try:
            import yaml

            cfg = yaml.safe_load(manifest.open())
            hits = [e["file"] for e in cfg["sources"]["ModelSEED"]["pkas"]
                    if e.get("consumed_by_production")]
            chosen = hits[0] if len(hits) == 1 else None
        except ImportError:
            # PyYAML is absent from some of the pinned environments. The
            # manifest is regular enough to scan without it, and doing so is
            # much better than silently reading a stale file: the fallback
            # export covers 22,399 compounds against the current run's 32,201.
            block, current = False, None
            for line in manifest.read_text().splitlines():
                if re.match(r"^  \w", line):
                    block = line.strip().startswith("ModelSEED:")
                if not block:
                    continue
                m = re.search(r"file:\s*(\S+)", line)
                if m:
                    current = m.group(1)
                if "consumed_by_production: true" in line and current:
                    chosen = current
        except Exception as exc:      # malformed manifest -- say so, do not guess
            raise RuntimeError(
                f"could not read {manifest}: {type(exc).__name__}: {exc}"
            ) from exc

    if chosen:
        path = base / Path(chosen).name
        if not path.exists():
            raise FileNotFoundError(
                f"{manifest} marks {chosen} as consumed_by_production but "
                f"{path} is missing; sync the snapshot from the repository"
            )
        return _read_pka_table(path)

    fallback = base / "opam2_molgpka.tsv"
    print(f"WARNING: no production pKa table resolved from {manifest}; "
          f"falling back to {fallback.name}, which is the superseded export "
          f"(22,399 compounds, hydrogen-indexed)", file=sys.stderr)
    return _read_pka_table(fallback)


def load_structures(dev_dir: Path = DEV_DIR) -> Dict[str, dict]:
    """Read cpd-keyed pH-7 structures from the aggregate file.

    Returns {cpd_id: {smiles, inchi_key, inchi, formula, charge, aliases}}.
    """
    out: Dict[str, dict] = {}
    path = Path(dev_dir) / "Unique_ModelSEED_Structures.txt"
    with path.open() as handle:
        next(handle, None)  # header
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            cpd, kind, aliases, formula, charge, structure = parts[:6]
            rec = out.setdefault(
                cpd,
                {
                    "smiles": None,
                    "inchi_key": None,
                    "inchi": None,
                    "formula": formula,
                    "charge": charge,
                    "aliases": [a for a in aliases.split(";") if a],
                },
            )
            if kind == "SMILE":
                rec["smiles"] = structure
            elif kind == "InChIKey":
                rec["inchi_key"] = structure
            elif kind == "InChI":
                rec["inchi"] = structure
    return out


# Alberty's package keys reactants by his own symbol and carries no structures,
# so the mapping to ModelSEED has to be by hand. Only entries verified against
# the derived ladder are listed; anything unmapped simply falls through to the
# next source.
ALBERTY_TO_SEED = {
    "atp": "cpd00002", "adp": "cpd00008", "amp": "cpd00018", "pi": "cpd00009",
    "ppi": "cpd00012", "coA": "cpd00010", "succinylcoA": "cpd00078",
    "acetate": "cpd00029", "ammonia": "cpd00013", "sulfate": "cpd00048",
    "sulfite": "cpd00081", "nitrate": "cpd00209", "nitrite": "cpd00075",
    "succinate": "cpd00036", "fumarate": "cpd00106", "malate": "cpd00130",
    "citrate": "cpd00137", "oxalate": "cpd00180", "pep": "cpd00061",
    "glucose6phos": "cpd00079", "fructose6phos": "cpd00072",
    "fructose16phos": "cpd00290", "glycerol3phos": "cpd00080",
    "ribose5phos": "cpd00101", "phosphoglycerate3": "cpd00169",
    "adenine": "cpd00128", "adenosine": "cpd00182", "imp": "cpd00114",
    "itp": "cpd00068", "idp": "cpd00090", "prpp": "cpd00103",
}

# kJ/mol per pKa unit at 298.15 K
_RT_LN10 = 8.314462618e-3 * 298.15 * math.log(10)


def load_alberty_pkas(package: Path = None) -> Dict[str, List[float]]:
    """Derive macroscopic pKa ladders from Alberty's BasicBiochemData package.

    The package stores, per species, ``{dfG, dfH, charge, nH}``. Consecutive
    species differ by one proton, so each adjacent pair gives a dissociation:
    ``pKa = [dfG(deprotonated) - dfG(protonated)] / (RT ln10)``. This is
    Alberty's own ``calcpK``.

    The values are at **zero ionic strength**, which is what eQuilibrator
    expects -- it applies Debye-Huckel itself from the microspecies charge --
    so no conversion is applied here. Validated against literature: ATP 7.60,
    Pi 7.22, PPi 9.46, acetate 4.75, ammonia 9.25, succinate 5.64/4.21.

    Ladders are truncated to the physiologically interesting range, because
    Alberty only tabulates species that matter between roughly pH 5 and 9;
    phosphate yields 7.22 alone, without its 2.15 and 12.35.
    """
    package = package or Path("/scratch/seaver/BasicBiochemData3.m")
    if not package.exists():
        return {}
    text = package.read_text(encoding="utf-8", errors="ignore")
    out: Dict[str, List[float]] = {}
    for name, body in re.findall(
        r"([A-Za-z][A-Za-z0-9]*)sp\s*=\s*(\{\{.*?\}\})\s*;", text, re.S
    ):
        cpd = ALBERTY_TO_SEED.get(name)
        if not cpd:
            continue
        species = []
        for row in re.findall(r"\{([^{}]*)\}", body):
            field = [x.strip() for x in row.split(",")]
            if len(field) != 4:
                continue
            try:                       # '_' marks an unmeasured enthalpy
                species.append((float(field[0]), int(float(field[3]))))
            except ValueError:
                continue
        species.sort(key=lambda t: t[1])
        ladder = [
            (a[0] - b[0]) / _RT_LN10
            for a, b in zip(species, species[1:])
            if b[1] - a[1] == 1
        ]
        if ladder:
            out[cpd] = sorted(ladder, reverse=True)
    return out


def load_iupac_pkas(table: Path = None) -> Dict[str, List[float]]:
    """Read the local IUPAC pKa flat file, if present.

    The file is produced by ``Scripts/Structures/Fetch_IUPAC_pKas.py`` and is
    **gitignored**: the dataset is CC-BY-NC-4.0 and is not redistributed with
    this database. Its absence is normal on a fresh clone and simply means the
    resolver falls through to MolGpKa.

    Values carry a bounded ionic-strength offset -- the dataset has no ionic
    strength column and most measurements were made in supporting electrolyte,
    so ATP reads 6.53 against Alberty's 7.60 at I=0. Measured against Alberty
    the mean error is 0.46 pKa units, against MolGpKa's 1.35, so it is used in
    preference to prediction but not in preference to Alberty.
    """
    if table is None:
        table = (Path(__file__).resolve().parents[3] / "Biochemistry"
                 / "Structures" / "ModelSEED" / "pkas" / "iupac_v2_3b.tsv")
    if not table.exists():
        return {}
    out: Dict[str, List[float]] = {}
    with table.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                out.setdefault(row["seed_id"], []).append(float(row["pka_value"]))
            except (TypeError, ValueError):
                continue
    return {k: sorted(v, reverse=True) for k, v in out.items()}


def load_cache_ladders(
    mapping: Path = None, cache: Path = None
) -> Dict[str, List[float]]:
    """Macroscopic ladders already held in the upstream cache, by ModelSEED id.

    ``build_seed_mapping.py`` has already matched our structures to cache rows
    on skeleton+stereo -- that is, ignoring only the protonation block -- and
    records the hit as ``our_cache_id``. This reads the ladder from each matched
    row, so no new matching is done here and the mapping stays the single place
    where structural identity is decided.

    Only classifications that assert the same molecule are honoured. A
    ``wrong_molecule`` or ``stereo_conflict`` row may still carry a cache id for
    reporting purposes, and taking its ladder would attach one compound's
    protonation to another -- the IPP/DMAPP failure mode.
    """
    mapping = mapping or (EQ_DIR / "data" / "seed_mapping.tsv")
    # MUST default to the pristine upstream cache. ``our_cache_id`` is an id in
    # THAT database (build_seed_mapping.py resolves it against
    # ~/.cache/equilibrator/compounds.sqlite), and reading a cache we ourselves
    # rebuilt would be circular -- it would hand back the MolGpKa values written
    # on the previous run and report them as an independent source.
    cache = cache or UPSTREAM_CACHE
    if not mapping.exists() or not cache.exists():
        return {}
    trusted = {"exact", "protonation_only", "remap", "regrouped_mobile_h"}
    wanted: Dict[int, List[str]] = {}
    with mapping.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = (row.get("our_cache_id") or "").strip()
            if cid and row.get("classification") in trusted:
                wanted.setdefault(int(cid), []).append(row["seed_id"])
    if not wanted:
        return {}
    con = sqlite3.connect(str(cache))
    out: Dict[str, List[float]] = {}
    for cid, seeds in wanted.items():
        row = con.execute(
            "select dissociation_constants from compounds where id=?", (cid,)
        ).fetchone()
        if not row or row[0] is None:
            continue
        v = row[0]
        try:
            v = pickle.loads(v) if isinstance(v, (bytes, bytearray)) else v
        except Exception:
            continue
        if not v:
            continue                   # empty is an answer, but not one to copy
        for seed in seeds:
            out[seed] = sorted((float(x) for x in v), reverse=True)
    con.close()
    return out


def resolve_pkas(
    structures: Dict[str, dict],
    marvin: Dict[str, Dict[str, List[float]]],
    molgpka: Dict[str, List[float]],
    preference: Tuple[str, ...] = DEFAULT_PKA_PREFERENCE,
    alberty: Dict[str, List[float]] = None,
    iupac: Dict[str, List[float]] = None,
    cache: Dict[str, List[float]] = None,
    exclusions: Optional[Set[str]] = None,
    admissible: Optional[Callable[[str, List[float]], bool]] = None,
) -> Dict[str, Tuple[List[float], str]]:
    """Pick one pKa list per compound and record where it came from.

    One source wins per compound; the sources are not merged. Mixing a
    macroscopic ladder with microscopic per-site values would produce a list
    that is neither, and eQuilibrator reads whatever it is given as a
    sequential ladder.

    An empty result is meaningful and distinct from absence: a compound the
    predictor ran on that reported no ionizable site has no pKas, and that is
    an answer. Callers should treat a missing key as "unknown" and an empty
    list as "does not ionize".

    Marvin values are keyed by source id, so they are reached through the
    compound's alias list; alberty and iupac are keyed by ModelSEED id
    directly, because both were matched on structure.
    """
    alberty = alberty or {}
    iupac = iupac or {}
    cache = cache or {}
    exclusions = load_substitution_exclusions() if exclusions is None else exclusions
    resolved: Dict[str, Tuple[List[float], str]] = {}
    for cpd, rec in structures.items():
        alias_set = set(rec["aliases"])
        if is_bare_ion(rec):
            # Assert the empty list rather than letting a fallback invent one.
            # NOTE: this stops us GENERATING the artefact; it does not clear a
            # value already inherited from the upstream cache, because the
            # builder treats an empty list as "no source of ours" and leaves the
            # existing row alone. Clearing those is task #49 and needs the
            # builder to distinguish empty-as-answer from absent.
            resolved[cpd] = ([], "none:monatomic")
            continue
        gated_ok = (cpd not in exclusions
                    and needs_macroscopic_fallback(cpd, molgpka))
        for tool in preference:
            if tool in GATED_TIERS and not gated_ok:
                continue        # MolGpKa handles this compound; do not displace it
            # A ladder that the magnesium guard will refuse is worse than no
            # ladder from this tier: the builder would skip re-protonation and
            # silently leave whatever was already stored. Falling through to the
            # next source instead is what lets phosphate keep a usable ladder
            # (Alberty's [7.22] alone spans too few proton counts) and lets PEP
            # take IUPAC's [6.38, 3.45] over Alberty's inadmissible [7.00].
            if tool == "alberty" and cpd in alberty:
                values = normalize_pkas(alberty[cpd])
                if values and (admissible is None or admissible(cpd, values)):
                    resolved[cpd] = (values, "alberty")
                    break
            elif tool == "iupac" and cpd in iupac:
                values = normalize_pkas(iupac[cpd])
                if values and (admissible is None or admissible(cpd, values)):
                    resolved[cpd] = (values, "iupac")
                    break
            elif tool == "cache" and cpd in cache:
                values = normalize_pkas(cache[cpd])
                if values and (admissible is None or admissible(cpd, values)):
                    resolved[cpd] = (values, "cache")
                    break
            if tool == "marvin":
                hit = None
                for source in MARVIN_SOURCES:
                    table = marvin.get(source, {})
                    for alias in alias_set:
                        if alias in table:
                            hit = (table[alias], f"marvin:{source}:{alias}")
                            break
                    if hit:
                        break
                if hit:
                    values = normalize_pkas(hit[0])
                    if values:
                        resolved[cpd] = (values, hit[1])
                        break
            elif tool == "molgpka" and cpd in molgpka:
                values = normalize_pkas(molgpka[cpd])
                if values:
                    # Last tier: accepted even if inadmissible, because the
                    # alternative is no ladder at all. The builder's guard still
                    # applies and will decline to install it.
                    resolved[cpd] = (values, "molgpka")
                    break
    return resolved


def magnesium_guard(cache_path: Path):
    """Build the ``admissible`` predicate ``resolve_pkas`` takes.

    Mirrors the guard in ``build_modelseed_cache.py``. Magnesium species are
    indexed by absolute proton count and were computed against the stored
    ladder, so a replacement ladder has to span the proton counts they
    reference or ``populate_microspecies`` fails looking for a reference that no
    longer exists.

    Reduced from the guard's set containment: with ``n_h`` the ``atom_bag``
    proton count and ``mg`` the referenced counts, a ladder is admissible iff it
    carries at least ``max(0, n_h - min(mg))`` values above pH 7 and at least
    ``max(0, max(mg) - n_h)`` at or below it. Compounds with no magnesium data
    are unconstrained.
    """
    facts: Dict[str, Tuple[int, List[int]]] = {}
    if cache_path and Path(cache_path).exists():
        con = sqlite3.connect(str(cache_path))
        seeds = cache_seed_identifiers(Path(cache_path))
        for seed, (cid, _ik) in seeds.items():
            row = con.execute(
                "select atom_bag from compounds where id=?", (cid,)).fetchone()
            if not row:
                continue
            bag = row[0]
            try:
                bag = pickle.loads(bag) if isinstance(bag, (bytes, bytearray)) else bag
            except Exception:
                bag = {}
            n_h = bag.get("H", 0) if isinstance(bag, dict) else 0
            mg = sorted({r[0] for r in con.execute(
                "select number_protons from magnesium_dissociation_constant "
                "where compound_id=?", (cid,))})
            if mg:
                facts[seed] = (n_h, mg)
        con.close()

    def admissible(cpd: str, values: List[float]) -> bool:
        got = facts.get(cpd)
        if not got:
            return True
        n_h, mg = got
        major = sum(1 for v in values if v > 7.0)
        reachable = {i - major + n_h for i in range(len(values) + 1)}
        return set(mg) <= reachable

    return admissible


def build_dataset(
    dev_dir: Path = DEV_DIR, preference: Tuple[str, ...] = DEFAULT_PKA_PREFERENCE
) -> Dict[str, dict]:
    """Join structures to pKas, keeping compounds that have both.

    Returns {cpd_id: {smiles, inchi_key, pkas, pka_source, ...}}.
    """
    structures = load_structures(dev_dir)
    resolved = resolve_pkas(
        structures, load_marvin_pkas(dev_dir), load_molgpka_pkas(dev_dir), preference
    )

    dataset = {}
    for cpd, rec in structures.items():
        if not rec["smiles"] or cpd not in resolved:
            continue
        pkas, source = resolved[cpd]
        entry = dict(rec)
        entry["pkas"] = pkas
        entry["pka_source"] = source
        dataset[cpd] = entry
    return dataset


_DECOMPOSER = None


def _get_decomposer():
    """Build the group decomposer once; it is expensive to construct."""
    global _DECOMPOSER
    if _DECOMPOSER is None:
        from equilibrator_assets import group_decompose

        _DECOMPOSER = group_decompose.GroupDecomposer()
    return _DECOMPOSER


def is_group_decomposable(smiles: str) -> bool:
    """Whether eQuilibrator's group decomposer can handle this structure."""
    from equilibrator_assets import molecule

    decomposer = _get_decomposer()
    try:
        mol = molecule.Molecule.FromSmiles(smiles)
        decomposer.Decompose(mol, ignore_protonations=False, raise_exception=True)
        return True
    except Exception:
        # GroupDecompositionError for unsupported chemistry, but OpenBabel can
        # also fail outright on a malformed structure.
        return False


def inchi_key_for(smiles: str) -> Optional[str]:
    """Compute an InChIKey the same way ``create_compounds`` does."""
    from openbabel.pybel import readstring

    try:
        return readstring("smiles", smiles).write("inchikey").strip()
    except Exception:
        return None


def pick_compound(candidates: Iterable[int], training: Set[int]) -> int:
    """Choose among cache compounds sharing a structure, as eQuilibrator does.

    ``get_or_create_compounds`` prefers a compound that is in the Component
    Contribution training set and otherwise falls back to the lowest id. This
    matters: 1,856 InChIKeys map to more than one compound in the shipped
    cache, and for 96 of them the training compound is *not* the highest id.
    Picking naively demotes common metabolites (water, phosphate, CO2, ammonia,
    pyruvate) from a measured reactant contribution to an unusable duplicate,
    which silently destroys every reaction that touches them.
    """
    candidates = sorted(candidates)
    for cid in candidates:
        if cid in training:
            return cid
    return candidates[0]


def cache_inchi_keys(
    sqlite_path: Path, training: Optional[Set[int]] = None
) -> Dict[str, int]:
    """Map every InChIKey in the compound cache to one internal compound id."""
    if training is None:
        training = training_compound_ids()
    con = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
    try:
        grouped: Dict[str, List[int]] = defaultdict(list)
        for cid, ik in con.execute(
            "select id, inchi_key from compounds where inchi_key is not null"
        ):
            grouped[ik].append(cid)
    finally:
        con.close()
    return {ik: pick_compound(cids, training) for ik, cids in grouped.items()}


def cache_seed_identifiers(sqlite_path: Path) -> Dict[str, Tuple[int, Optional[str]]]:
    """Current seed: accession -> (internal compound id, that compound's key).

    These are the MetaNetX-derived mappings that Path A replaces.
    """
    con = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
    try:
        return {
            acc: (cid, ik)
            for acc, cid, ik in con.execute(
                "select ci.accession, cp.id, cp.inchi_key "
                "from compound_identifiers ci "
                "join registries r on ci.registry_id = r.id "
                "join compounds cp on cp.id = ci.compound_id "
                "where r.namespace = 'seed'"
            )
        }
    finally:
        con.close()


def normalize_formula(formula: Optional[str]) -> Optional[str]:
    """Canonicalise a formula string for comparison (element-sorted, explicit 1s)."""
    import re

    if not formula or formula.strip() in ("null", "None", ""):
        return None
    parts = [
        (m[0], m[1] or "1")
        for m in re.findall(r"([A-Z][a-z]?)(\d*)", formula)
        if m[0]
    ]
    return "".join(f"{e}{n}" for e, n in sorted(parts)) or None


def cache_formulas(sqlite_path: Path) -> Dict[int, str]:
    """Molecular formula per cache compound, derived from its atom_bag.

    Used to separate genuine mapping errors (different formula) from
    representational differences such as ring forms and tautomers, which share
    a formula and should not be "corrected".
    """
    import pickle

    con = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
    out: Dict[int, str] = {}
    try:
        for cid, blob in con.execute(
            "select id, atom_bag from compounds where atom_bag is not null"
        ):
            try:
                bag = pickle.loads(blob)
            except Exception:
                continue
            if not bag:
                continue
            raw = "".join(f"{e}{n}" for e, n in sorted(bag.items()) if e != "e-")
            norm = normalize_formula(raw)
            if norm:
                out[cid] = norm
    finally:
        con.close()
    return out


def training_compound_ids() -> Set[int]:
    """The reactant-contribution compound ids baked into cc_params.npz."""
    from equilibrator_api import ComponentContribution

    return set(int(i) for i in ComponentContribution().predictor.preprocess._compound_ids)
