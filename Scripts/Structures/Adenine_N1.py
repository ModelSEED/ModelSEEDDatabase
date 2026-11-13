#!/usr/bin/env python
"""Adenine N1: Marvin 26.1's basic pKa there is wrong; use Marvin 23.4's.

THE DEFECT. Marvin 26.1 predicts the basic pKa of adenine N1 at 7.0-7.3
(adenosine 7.02, ATP 7.27, NAD 7.06, CoA 7.24 in pkas/marvin_26.1.tsv). At pH 7
that straddles the protonation threshold, so the 26.1 major microspecies
protonates N1 in 2,560 of the 2,914 adenine-containing SMILE rows of the
protonation bundles. Marvin 23.4 put the same site at 2.2-2.5 (adenosine 2.52,
ATP 2.21, NAD 2.50, CoA 2.53) and protonated N1 in 4 of 2,951 rows. The
measured value is about 3.5-4.2, so 23.4 is right and 26.1 is the outlier.
Nothing upstream of Marvin can fix this: RDKit already hands Marvin the
standard 9H-amino tautomer, and `cxcalc pka` on that tautomer still returns
~7.0 for N1.

THE OVERRIDE. 23.4's value is used in both places 26.1 is consumed:

  * protonation structures (Run_Marvin_Protonations.py): after the
    microspecies step, an adenine N1 carrying H and +1 is returned to neutral
    with one H fewer -- what 23.4's pKa (~2.4, far below 7) implies at pH 7.
    dH == dcharge == -1, so the invariant gate still holds;
  * pKa ladders (pkas/marvin_<ver>.tsv, Override_Adenine_N1_pKas.py): the
    pKb token on N1 is replaced by 23.4's N1 value for the same compound when
    23.4's ladder identifies it unambiguously (exactly one basic value in
    N1_234_RANGE), else by the median of those unambiguous values.

Only N1 is touched, and only where 26.1 has it protonated (structures) or
at or above N1_26_THRESHOLD (ladders). 1-methyladenine-type N1 (no H; a
permanent cation) is left alone, as are N3, N7 and the exocyclic N6.

ADENINE is matched as a 6-aminopurine with an acyclic, single-bonded N6, so
N6-substituted adenines (N6-methyl, N6-isopentenyl, zeatin) are covered and
2-amino/6-oxo purines (guanine, hypoxanthine) are not.
"""

import statistics

from rdkit import Chem

# N6 (acyclic, single bond) - C6 : N1 : C2 : N3 : C4 : C5 : N7 : C8 : N9.
# Atom 2 of each match is N1. Bond orders are left open (~) so that Kekule,
# aromatic, and charged forms written by Marvin all match.
ADENINE = Chem.MolFromSmarts(
    "[#7;!R]-[#6;R]1~[#7;R]~[#6;R]~[#7;R]~[#6;R]2~[#6;R]~1~[#7;R]~[#6;R]~[#7;R]~2")
N1 = 2

# 23.4's N1 values sit at 2.2-2.6; any other basic site in an adenine-containing
# compound that falls in this window makes the 23.4 ladder ambiguous.
N1_234_RANGE = (0.5, 5.0)
# A 26.1 N1 value below this is already on the right side of pH 7 and is kept.
N1_26_THRESHOLD = 5.0


def _rdkit_view(mol, MolExporter):
    """RDKit Mol over a Marvin Molecule, same atom order, unsanitized."""
    block = str(MolExporter.exportToFormat(mol, "mol"))
    rd = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)
    if rd is None:
        return None
    rd.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(rd)
    return rd


def n1_atoms(mol, MolExporter):
    """0-based indices of every adenine N1 in a Marvin Molecule."""
    rd = _rdkit_view(mol, MolExporter)
    if rd is None:
        return []
    return sorted({m[N1] for m in rd.GetSubstructMatches(ADENINE)})


def deprotonate_n1(mol, MolExporter):
    """Return adenine N1-H(+) to neutral N1 in place. Returns sites changed.

    Only an N1 that carries +1 AND at least one H is changed; an N1 bonded to
    carbon on all three sides (1-methyladenine) has no H and is a permanent
    cation, not a protonation state.
    """
    rd = _rdkit_view(mol, MolExporter)
    if rd is None:
        return 0
    targets = set()
    for match in rd.GetSubstructMatches(ADENINE):
        a = rd.GetAtomWithIdx(match[N1])
        if a.GetFormalCharge() == 1 and a.GetTotalNumHs(includeNeighbors=True) >= 1:
            targets.add(match[N1])
    # Highest index first, so removing an explicit H atom (always appended
    # after heavy atoms by Marvin) cannot shift a pending target's index.
    for i in sorted(targets, reverse=True):
        atom = mol.getAtom(i)
        explicit_h = [atom.getLigand(j) for j in range(atom.getBondCount())
                      if atom.getLigand(j).getAtno() == 1]
        if explicit_h:
            mol.removeAtom(explicit_h[0])
        else:
            # setCharge alone keeps the implicit H, leaving a neutral N1-H
            # (the imino tautomer) with one H too many.
            atom.setImplicitHcount(max(atom.getImplicitHcount() - 1, 0))
        atom.setCharge(0)
    return len(targets)


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


def n1_value_234(pkb_value):
    """23.4's N1 pKb from its pKb ladder string, or None if ambiguous/absent."""
    lo, hi = N1_234_RANGE
    hits = [v for _, _, v in parse_tokens(pkb_value) if lo <= v <= hi]
    return hits[0] if len(hits) == 1 else None


def median_n1_234(pkb_by_id):
    """Median of the unambiguous 23.4 N1 values over {id: pKb ladder string}."""
    vals = [v for v in (n1_value_234(s) for s in pkb_by_id.values()) if v is not None]
    return round(statistics.median(vals), 2) if vals else None
