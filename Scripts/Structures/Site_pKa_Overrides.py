#!/usr/bin/env python
"""Site classes where Marvin 26.1's pKa is measurably wrong, and their correction.

Two classes, both measured, both narrow:

  adenine_n1     6-aminopurine N1. 26.1 predicts the basic pKa at 7.0-7.3
                 (adenosine 7.02, ATP 7.27, NAD 7.06, CoA 7.24); Alberty's
                 tables (Structures/ModelSEED/pkas/literature_alberty-2003.tsv,
                 I = 0) give 4.68 ATP, 4.36 ADP, 3.99 AMP; Marvin 23.4 gave
                 2.2-2.5. Target: 2.53, the median of 23.4's unambiguous N1
                 values (at pH 7 any value below ~5 is equivalent).
  his_imidazole  4(5)-substituted imidazole on a CH2-C-N(amine) arm: the
                 histidine / histamine family, where the nearby ammonium
                 lowers the ring pKa. 26.1 predicts 7.20 (histidine) and 7.30
                 (histamine) against IUPAC-measured 6.04 and 6.05 (23.4: 6.14,
                 6.24). Free imidazole (26.1 7.09, measured 6.95), acylated
                 arms such as carnosine (26.1 6.94, measured ~6.8) and
                 N-substituted or metal-bound rings are NOT in the class: 26.1
                 is right there, or the site is not a protonation site.
                 Target: 6.1.

HOW THE STRUCTURES ARE CORRECTED (Run_Marvin_Protonations.py). Marvin's
microspecies distribution is reweighted: every microstate in which a class
site is protonated is multiplied by 10**(target - TYPICAL_261), the shift of
that site's microscopic pKa. The weights are then summed per net charge (the
macrospecies -- what a pKa ladder and eQuilibrator's major-species index
describe), the most populated charge is chosen, and within it the most
populated microstate. This MOVES a proton rather than dropping it: where
Marvin had parked the -3 state's proton on adenine N1 of ATP, the reweighted
distribution puts it on the terminal phosphate, as 26.1's own phosphate pKa
(7.57) and Alberty (7.60, I = 0) say. The earlier post-hoc N1 deprotonation
removed it and made ATP -4.

Protonation is judged by H count AND charge, never by one alone: a ring N
bonded to cobalt carries +1 and no proton, and the 2-thione tautomer of
ergothioneine has both ring N protonated at zero charge.

HOW THE LADDERS ARE CORRECTED (Override_Site_pKas.py). The 26.1 pKb token on
a class site is replaced by 23.4's value for the same compound where 23.4's
ladder identifies it unambiguously (exactly one basic value in the class
window), else by the median of those values over all four sources.

The reweighting uses class constants rather than per-compound 26.1 values:
both classes are narrow enough that the 26.1 prediction is tight (N1
7.0-7.3; His/histamine 7.20-7.30), and for N1 any shift below ~-3 already
makes the protonated microstates negligible.
"""

import statistics

from rdkit import Chem

METALS = {"Co", "Fe", "Ni", "Mg", "Zn", "Cu", "Mn", "Mo", "W", "Ca", "Cd", "Hg", "V", "Cr"}

CLASSES = {
    "adenine_n1": dict(
        # N6 (acyclic, single bond) - C6 : N1 : C2 : N3 : C4 : C5 : N7 : C8 : N9
        smarts="[#7;!R]-[#6;R]1~[#7;R]~[#6;R]~[#7;R]~[#6;R]2~[#6;R]~1~[#7;R]~[#6;R]~[#7;R]~2",
        sites=(2,),             # N1
        target=2.53, typical_261=7.1,
        window_234=(0.5, 5.0),  # 23.4's N1 value is the single basic value in here
        threshold_261=5.0,      # a 26.1 ladder token at or above this is replaced
    ),
    "metal_macrocycle": dict(
        # Corrin / porphyrin / chlorophyll: a metal (Co/Mg/Fe/...) anywhere in the
        # molecule plus a macrocyclic ring (>=12 atoms, >=3 N). Sites are the
        # ring N atoms; the protonated one carries the H/+1. Target taken from
        # 23.4's ladder via window_234, NOT a class constant, because corrin vs
        # porphyrin N differ.
        smarts=None,            # matched programmatically in find_sites (metal + big ring)
        sites=None,
        target=None, typical_261=7.6,
        window_234=(5.0, 8.5),
        threshold_261=6.9,
        macrocycle=True,
    ),
    "heterocycle_n": dict(
        # Any aromatic ring N (5- or 6-membered) not covered by the two classes
        # above: purines other than adenine N1, pteridines, benzimidazoles,
        # drug heterocycles. Sites = all aromatic ring N atoms.
        smarts="[n,N;R]",              # any ring N, aromatic or aliphatic (covers alkaloid N+ too)
        sites=None,             # all matches
        target=None, typical_261=7.2,
        window_234=(4.0, 6.9),
        threshold_261=6.9,
    ),
    "his_imidazole": dict(
        # ring C(arm) : C : N : C2 : N, arm CH2-C-N with N not acylated/sulfonylated.
        # C2 must be an unsubstituted CH (D2): 2-thio/2-oxo/2-amino imidazoles
        # (ergothioneine, thiolhistidine) are thiones/ureas, not imidazoles.
        smarts=("[#6;R1;r5]1(-[CH2]-[CX4]-[#7;!$([#7]-[#6]=[#8]);!$([#7]-[#16])])"
                "~[#6;R1;r5]~[#7;R1;r5;D2]~[#6;R1;r5;D2]~[#7;R1;r5;D2]~1"),
        # match order: 0 ring C, 1 CH2, 2 C, 3 arm N, 4 ring C, 5 ring N, 6 ring C, 7 ring N
        sites=(5, 7),           # the two ring N (both unsubstituted: D2)
        target=6.1, typical_261=7.25,
        window_234=(5.0, 6.9),
        threshold_261=6.9,
    ),
}
_PATTERNS = {}


def pattern(name):
    if name not in _PATTERNS:
        _PATTERNS[name] = (Chem.MolFromSmarts(CLASSES[name]["smarts"])
                           if CLASSES[name]["smarts"] else None)
    return _PATTERNS[name]


def macrocycle_sites(rd):
    """Ring N atoms of a metal-bound macrocycle (>=12-membered ring with >=3 N,
    metal anywhere in the molecule). Corrin/porphyrin/chlorophyll central N."""
    out = []
    for ri in rd.GetRingInfo().AtomRings():
        if len(ri) >= 12:
            ns = tuple(sorted(i for i in ri if rd.GetAtomWithIdx(i).GetSymbol() == "N"))
            if len(ns) >= 3:
                out.append(ns)
    return out


def rdkit_view(mol, MolExporter):
    """RDKit Mol over a Marvin Molecule, same atom order, unsanitized."""
    rd = Chem.MolFromMolBlock(str(MolExporter.exportToFormat(mol, "mol")),
                              sanitize=False, removeHs=False)
    if rd is None:
        return None
    rd.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(rd)
    return rd


def find_sites(rd):
    """[(class, (atom indices of the site))] in an RDKit view; metal fragments skipped
    for the imidazole class (a ring N bound to a metal is not a protonation site)."""
    if rd is None:
        return []
    has_metal = any(a.GetSymbol() in METALS for a in rd.GetAtoms())
    out = []
    for name, spec in CLASSES.items():
        if name == "his_imidazole" and has_metal:
            continue
        if name == "metal_macrocycle":
            if not has_metal:
                continue
            for site in macrocycle_sites(rd):
                out.append((name, site))
            continue
        seen = set()
        pat = pattern(name)
        matches = rd.GetSubstructMatches(pat) if pat else ()
        for m in matches:
            if spec["sites"]:
                site = tuple(sorted(m[i] for i in spec["sites"]))
            else:
                site = tuple(sorted(m))
            if site not in seen:
                seen.add(site)
                out.append((name, site))
    return out


def _h(atom):
    return atom.getImplicitHcount() + atom.getExplicitHcount()


def site_protonated(marvin_mol, name, site):
    """Is this class site protonated in a Marvin microspecies? By H count."""
    atoms = [marvin_mol.getAtom(i) for i in site]
    if name == "adenine_n1":
        a = atoms[0]
        return a.getCharge() == 1 and _h(a) >= 1
    # imidazolium: both ring N carry H AND the ring N carry +1 in total.
    # H count alone would also catch a 2-thione tautomer (both N-H, neutral);
    # charge alone would catch a metal-bound [n+] with no proton.
    return (sum(_h(a) for a in atoms) >= 2
            and sum(a.getCharge() for a in atoms) == 1)


def reweighted_major(plugin, sites, site_targets=None):
    """Major microspecies after the class reweighting. Assumes plugin.run() was
    called on the fragment the sites were found in (microspecies keep its atom
    order). Returns (molecule, {class: bool changed vs Marvin's own pick}).

    site_targets maps (class, site) -> per-compound target pKa for classes whose
    CLASSES target is None (metal macrocycles, generic heterocycles); it is
    derived from 23.4's ladder by the caller."""
    site_targets = site_targets or {}
    n = plugin.getMicrospeciesCount()
    rows = []
    for i in range(n):
        m = plugin.getSortedMicrospecies(i)
        w = plugin.getSortedMsDistribution(i)
        f = 1.0
        hit = set()
        for name, site in sites:
            if site_protonated(m, name, site):
                spec = CLASSES[name]
                target = spec["target"]
                if target is None:
                    # per-compound target from 23.4's ladder (set by the caller
                    # via site_targets); fall back to typical_261 shift of -1
                    target = site_targets.get((name, site), spec["typical_261"] - 1.0)
                f *= 10 ** (target - spec["typical_261"])
                hit.add(name)
        q = sum(m.getAtom(j).getCharge() for j in range(m.getAtomCount()))
        rows.append((w * f, q, m, hit))
    if not rows:
        return plugin.getMajorMicrospecies(), {}
    by_charge = {}
    for wf, q, _, _ in rows:
        by_charge[q] = by_charge.get(q, 0.0) + wf
    best_q = max(by_charge, key=by_charge.get)
    best = max((r for r in rows if r[1] == best_q), key=lambda r: r[0])
    marvin_pick = rows[0]                         # sorted: index 0 is Marvin's major
    changed = {name: (marvin_pick[2] is not best[2]) and (name in marvin_pick[3] or name in best[3])
               for name, _ in sites}
    return best[2], changed


# ------------------------------------------------------------------ ladders
def parse_tokens(value):
    """'f:a:v;f:a:v' -> [(fragment, atom, value)] with numeric value."""
    out = []
    for tok in (value or "").split(";"):
        parts = tok.split(":")
        if len(parts) == 3:
            out.append((parts[0], parts[1], float(parts[2])))
    return out


def format_tokens(tokens):
    return ";".join(f"{f}:{a}:{v:.2f}" for f, a, v in tokens)


def value_234(name, pkb_value):
    """23.4's value for this class from its pKb ladder string, or None if ambiguous/absent."""
    lo, hi = CLASSES[name]["window_234"]
    hits = [v for _, _, v in parse_tokens(pkb_value) if lo <= v <= hi]
    return hits[0] if len(hits) == 1 else None


def median_234(name, pkb_strings):
    vals = [v for v in (value_234(name, s) for s in pkb_strings) if v is not None]
    if vals:
        return round(statistics.median(vals), 2)
    t = CLASSES[name]["target"]
    return t if t is not None else 6.0   # generic fallback for target=None classes
