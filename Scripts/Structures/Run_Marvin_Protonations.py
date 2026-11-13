#!/usr/bin/env python
"""Regenerate the per-source Marvin protonation bundles at pH 7.

Replaces the ChargeMol.java path (which needed MarvinBeans on the classpath)
with chemaxon.calculations.MajorMicrospeciesPlugin driven through JPype.

  input   Biochemistry/Structures/<source>/inchi.tsv    (preferred)
          Biochemistry/Structures/<source>/smiles.tsv   (gap-fill only)
  output  Biochemistry/Structures/<source>/protonations/marvin_<ver>_ph<n>.tsv
          columns: external_id, type, structure, formula, charge,
                   tool, tool_version, ph, generated_on

METHODOLOGY DELTA -- READ THIS FIRST. ChargeMol.java, which produced the 23.4
bundles, ran FIVE steps per fragment:

    i)   convert explicit hydrogens to implicit
    ii)  aromatize
    iii) find the dominant tautomer at pH 7
    iv)  find the major microspecies at pH 7
    v)   convert implicit hydrogens back to explicit

Step (iii) was originally chemaxon TautomerizationPlugin (Isomers Plugin
Group). The licence that produced the first 26.1 bundles covered only the
Protonation Plugin Group (pKa + major microspecies); Isomers was refused:

    $ cxcalc "CC(=O)CC(=O)C" majortautomer -H 7
    chemaxon.license.api.LicenseException: No valid license has been found.
    Product name: Isomers Plugin Group

So those bundles ran (i),(ii),(iv),(v) and SKIPPED (iii) — see the archived
Marvin_26.1_Protonation_Run.md for the 20.3% net-charge impact.  This script
now restores (iii) via RDKit's canonical tautomer
(rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumerator,
Sitzmann et al. 2010, SetRemoveSp3Stereo/SetRemoveBondStereo=False to
preserve chiral centres adjacent to tautomerizable bonds) — see the
module-level tautomer_canonical() for the full contract.  Two review
fixes are baked in: the input is first neutralized with
rdMolStandardize.Uncharger, because MetaCyc/ChEBI/Rhea store charged
species whose InChIs RDKit parses as imidate ANIONS N=C([O-]), and
tautomer enumeration never moves charge, so without uncharging the
imidate defect survives; this is safe because the Marvin microspecies
step recomputes the pH-7 charge and the formula/charge gate still checks
the final row against the source.  And when the enumerator hits its
MaxTautomers/MaxTransforms cap (the only wall-clock bound that can
interrupt the C++ enumeration) the partial canonical pick is STILL used
and the cap counted as tautomer_cap_hit, rather than abandoning the
compound.  RDKit's scoring is a fixed aromaticity/carbonyl rule, not
pKa-aware, so this is NOT a pH-7 dominant-tautomer in the Marvin sense;
parity against Marvin 23.4's shipped output is measured by
Scripts/Structures/Benchmark_RDKit_Tautomer.py (report:
Biochemistry/Structures/_reports/benchmark_rdkit_tautomer.tsv) on the
compounds where 26.1-without-(iii) disagreed with 23.4.  The tautomer
step is now (iii-RDKit) between aromatize and microspecies; nothing else
changes.  Bundle provenance records both tool versions (Marvin 26.1 +
RDKit) and tautomer method.

ADENINE N1 is the one site where the 26.1 engine itself, not the missing
tautomer step, disagrees with 23.4. 26.1 puts its basic pKa at 7.0-7.3
(23.4: 2.2-2.5; measured ~3.5-4.2), so the microspecies protonates N1 in
2,560 of 2,914 adenine rows (23.4: 4 of 2,951) and every adenine nucleotide,
NAD(P), FAD and CoA ester is one unit too positive. RDKit already hands Marvin
the standard amino tautomer, so step (iii) cannot help. After the microspecies
step an adenine N1-H(+) is returned to neutral N1 -- 23.4's value applied at
pH 7 -- and counted as adenine_n1_override. dH == dcharge == -1, so the
invariant gate is unaffected. The pKa ladders get the same correction from
Override_Adenine_N1_pKas.py. See Adenine_N1.py. Measured end-to-end through
Marvin 26.1 on the 3,000-compound benchmark sample: charge parity with 23.4
51.9% -> 72.8%; all 394 CoA thioesters match 23.4.

INPUT SELECTION follows Run_Marvin_pKas.py: `smiles.tsv` is a strict superset
of `inchi.tsv`, and InChI wins wherever it exists, because a charged SMILES is
an already-deprotonated species and Marvin would be answering a different
question than it does for the neutral InChI parent. Each compound is protonated
exactly once, from its preferred representation.

WITH ONE EXCEPTION, and it is not a small one. InChI DISCONNECTS METAL-LIGAND
BONDS by design. Triphenyltin chloride is stored in smiles.tsv as the intact
molecule `Cl[Sn](c1ccccc1)(c1ccccc1)c1ccccc1`, but its InChI is

    InChI=1S/3C6H5.ClH.Sn/c3*1-2-4-6-5-3-1;;/h3*1-5H;1H;/q;;;;+1/p-1

-- three phenyl radicals, HCl and a tin atom, as five separate components. Take
InChI-first on that and Marvin is handed an already-shattered molecule; it
protonates each piece on its own and the bundle ships
`[Cl-].[SnH3+].[c]1ccccc1.[c]1ccccc1.[c]1ccccc1` where 23.4 shipped the intact
structure. 487 compounds across the four sources have an InChI more fragmented
than their SMILES -- cobalamins, Ni/Fe/Mg porphyrins, molybdenum cofactors,
organotins -- and the first cut of this bundle shattered 461 of them.

So the rule is InChI-first UNLESS the InChI is the more fragmented of the two,
in which case the SMILES wins. Comparing fragment counts rather than screening
for metals keeps this general: whatever the reason an InChI has taken a
molecule apart, the representation that keeps it together is the better input.
Counted as `smiles_preferred` in the per-source stats.

For those same compounds the InChI ROW is built from the InChI-derived
molecule rather than the SMILES-derived one (`inchi_from_inchi_form`), so each
column carries what its own representation can express. It is PROTONATED, not
passed through: 23.4 protonated this form uniformly, and the protonation merely
happened to be a no-op for 234 of the 493 while changing the other 259 -- the
chlorophylls and cobalamins pick up a /p-2 layer. Passing the source through
would silently un-protonate those 259.

The cost is visible and is recorded rather than hidden: Marvin 26.1 reads the
detached `4Fe.4S` as four FREE sulfide ions and protonates them to H2S at pH 7,
where 23.4 left them alone, so ChEBI 33722 ships Fe4S4 as H8Fe4S4 and 136511
ships MnO2 as H4MnO2. That is the same class of engine change as the 20.3%
net-charge delta above, confined to ligands InChI detached from their metal.

Nothing in this is RDKit's doing, and swapping the importers would not help:
Marvin reads the same disconnected InChI the same way. The loss happens in the
InChI string itself, before any parser sees it.

OUTPUT LAYOUT reproduces the 23.4 bundles exactly, which is one row per
representation rather than one row per compound:

  SMILE             every compound in smiles.tsv
  InChI, InChIKey   only compounds that also appear in inchi.tsv

The asymmetry is chemistry, not omission. The compounds with no InChI are the
8,834 polymers (`*` attachment points) and organometallics; InChI cannot
represent a query molecule, so those compounds carry a SMILE row and nothing
else -- exactly as they did under 23.4, where the SMILE-only ids match the
structure-file smiles-only ids one for one in every source.

ONE ENGINE. Every compound goes through the Java plugin, including the ones the
cxcalc CLI could have handled. This is deliberate: the three representations
must describe the SAME protonated molecule, and exporting all three from one
in-memory Molecule guarantees that.

The plugin agrees with the CLI on the chemistry: across 300 ChEBI compounds,
formula and charge are identical in 296. All 4 exceptions are organometallics
(Mg/Fe/Ni porphyrins and chlorophylls) where `cxcalc majorms` silently DROPS
the metal fragment -- it returns C35H34N4O5 for a compound whose protonated
form is C35H34MgN4O5 -- and the per-fragment route here keeps it. That is the
same reason ChargeMol.java fragmented and re-fused rather than protonating the
molecule whole, and it is why the CLI is not used for the bulk.

TWO LOSSY IMPORTERS, AND THE BRIDGE AROUND THEM. Marvin's InChI and SMILES
importers do not preserve stereochemistry. Measured by round-tripping the
source InChI with no protonation at all -- import, export, compare:

    MolImporter.importMol(<InChI>)        ->  61% reproduce the source
    MolImporter.importMol(<SMILES>)       ->  39%
    MolImporter.importMol(<RDKit molblock>) -> 99.6%
    RDKit alone, InChI -> mol -> InChI     -> 99.7%   (control)

The loss is on IMPORT, not on export: both Marvin writers score the same ~60%
behind a lossy import, and the RDKit control proves the InChI itself carries
everything needed. So every structure is bridged through RDKit into an MDL
molblock before Marvin sees it, which is what ChargeMol.java was doing when it
read `args[0]` as a mol file rather than taking a string. 2D coordinates are
not required -- the molblock's parity flags are enough.

This matters against the benchmark: where 23.4 protonation changed nothing,
23.4's InChI reproduces the source in 96.6% of ChEBI compounds. An unbridged
run scores 55-61% and would have shipped a stereochemistry regression in tens
of thousands of rows.

EXPORT uses `inchi:AuxNone,SAbs`. The SAbs flag (absolute stereo) is not
cosmetic: without it the same bridged molecule reproduces the source InChI only
30.9% of the time, against 99.6% with it.

The InChIKey is DERIVED from that InChI string rather than exported separately,
because Marvin's `inchikey` export does not honour SAbs and so disagrees with
its own InChI export. For Rhea POLYMER_10033 the two exports describe the same
molecule but hash differently -- Marvin's key is ...-ARPYZQPTNA-N while the key
of the InChI it just wrote is ...-VFUOTHLCSA-N, which is what 23.4 shipped.
Deriving makes the two columns consistent by construction.

The plugin does place charge on wildcard-adjacent atoms -- on the Rhea glycans,
the peptide-backbone nitrogen hidden behind `*` is protonated to [NH2+] because
Marvin cannot see the amide bond the wildcard stands for. This is inherited,
not introduced: 23.4 wrote `*[NH2+][C@@H](...` for POLYMER_12621 and this run
writes the same. Correcting it would be a curation decision about what `*`
means, not a regeneration, so it is left alone and recorded here.

FAILURE TAXONOMY. A compound yields no row only when no rung of the fallback
ladder in protonate_best() produces a writable molecule. Across all four
sources that leaves two categories, and only one of them is a failure:

  H+ itself (`InChI=1S/p+1`, ChEBI 15378 / KEGG C00080) has no heavy atom and
  so has no protonated form. 23.4 carries no row for it either; agreeing with
  the old bundle here is the correct outcome, not a loss.

  ChEBI 60492 is the dative-bond SMILES (`...[Mg]35<-N2=...`) whose `<`
  character the SMILES parser rejects. This is the same single structure the
  pKa regeneration could not read, and it is a genuine loss of one row against
  23.4. Recording the id rather than only the count, because the last
  regeneration was reviewed for reporting failures as a bare number.

Everything else that the first attempt could not write -- query molecules
broken by aromatize(), and InChIs like azide's that rebuild as radicals -- is
recovered by the ladder rather than dropped. See protonate_best().

FORMULA AND CHARGE come from Print_Structure_Formula_Charge.parse_structure --
this repository's own function, imported rather than reimplemented, computed per
ROW from that row's structure string.

Do NOT substitute Marvin's getFormula() here. Marvin omits wildcard atoms from
a formula; this repository renders them as R, in one line of that function:

    formula = re.sub(r'\\*', 'R', formula)

The first cut of this script wrote Marvin's formula and deferred the refresh to
Print_Structure_Formula_Charge.py as a follow-up step. That was wrong, and it
shipped: every one of the 8,704 rows that should carry an R group lost it,
taking the count of SMILE rows whose formula contains R from 23.4's 8,704 down
to zero. Stearoyl-ACPs went from
C32H60N3O9PR2S to C32H60N3O9PS -- and downstream, `Update_Compound_Structures_
Formulas_Charge.py` propagated that into 6,052 compound records, turning
cpd00049 "carboxylic acid" from CHO2R into CHO2. A generic compound stopped
being generic.

The STRUCTURES were never affected -- 8,727 SMILE structures carry a `*` in
both bundles, identically -- which is exactly why this was invisible in every
coverage count and every InChI comparison. Only the formula column was wrong,
and the cascade consumes that column directly, so the bundle has to be correct
as written rather than correct after a second script runs.

Marvin's getFormula() survives as the fallback for the handful of rows
parse_structure cannot read, counted as `formula_from_marvin` in the per-source
stats so it can never again be a silent substitution.

When checking any of this, count R with HAS_R_GROUP (`R(?![a-z])`) and not a
substring test for "R": Ru, Rb, Rh, Re and Rn all match the naive test and
inflate the count by 8 across these bundles. The first report of this
regression said 8,735 and 8,712 for that reason, which made it look as though
eight rows had survived when in fact none had.

REQUIRES. `pip install jpype1` and rdkit (already a dependency of
Print_Structure_Formula_Charge.py), plus the Marvin jars -- found next to the
cxcalc binary, or point MARVIN_LIB at <marvinsuite>/lib. A JRE is sufficient;
nothing here is compiled.

This script does NOT write compound records. Update_Compound_Structures_
Formulas_Charge.py is the step that rewrites them.
"""

SOURCES = ["ChEBI", "KEGG", "MetaCyc", "Rhea"]

if __name__ == "__main__":
    # Argument guard -- see "The argument guard" in Scripts/README.md.
    import argparse as _argparse
    _parser = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter)
    # No choices= here: with nargs="*" argparse validates the DEFAULT against
    # choices as a single value, so any default fails. Validated by hand below.
    _parser.add_argument(
        "sources", nargs="*", default=[],
        help=f"sources to regenerate, from {SOURCES} (default: all four)")
    _parser.add_argument(
        "--ph", type=float, default=7.0,
        help="pH for the major-microspecies calculation (default: 7.0)")
    _parser.add_argument(
        "--limit", type=int, default=0,
        help="process at most N compounds per source (smoke tests)")
    _ARGS = _parser.parse_args()
    _bad = [s for s in _ARGS.sources if s not in SOURCES]
    if _bad:
        _parser.error(f"invalid source(s) {_bad}; choose from {SOURCES}")


import glob
import os
import re
import shutil
import subprocess
import sys
from csv import DictReader
from datetime import date

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "Libs", "Python"))
# The formula/charge columns of this file belong to Print_Structure_Formula_Charge,
# so they are computed with ITS function rather than a second implementation --
# see FORMULA AND CHARGE.
from Print_Structure_Formula_Charge import parse_structure    # noqa: E402
from BiochemPy import Compounds                              # noqa: E402

_RDKIT_VER_STR = None

def _rdk_ver_str(version):
    """tool_version string with RDKit provenance, computed once."""
    global _RDKIT_VER_STR
    if _RDKIT_VER_STR is None:
        try:
            from rdkit import __version__ as _v
        except Exception:
            _v = ""
        _RDKIT_VER_STR = f"{version}+rdkit{_v}" if _v else version
    return _RDKIT_VER_STR
# THE GATE. A protonation moves protons, so a protonated row must satisfy
# dH == dcharge against its source with heavy atoms conserved. A row that does
# not is not a protonation state -- InChI disconnected a metal cluster and the
# freed ligands were protonated as ions, or a bare S/Se/P/O atom came back as
# its hydride -- and it is replaced by its own source row. Same function the
# post-hoc Repair_Protonation_Rows.py applies, so a bundle repaired after the
# fact and one gated at generation are identical by construction.
import Validate_Protonations as _inv                         # noqa: E402
# Adenine N1 override (26.1 pKb ~7.0 vs 23.4 ~2.4) -- see ADENINE N1.
import Adenine_N1                                            # noqa: E402

# What Marvin calls a wildcard atom, across the notations it accepts. `R#` is
# the one that matters and the one easiest to miss: a molblock R atom -- which
# is how RDKit writes every dummy atom, and therefore how EVERY structure
# arrives here through the import bridge -- reads back as symbol "R#", not "R".
# A set without it counts zero wildcards on a molecule that plainly has them.
WILDCARD_SYMBOLS = {"*", "A", "R", "R#"}

# An R group in a formula, not the R of Ru/Rb/Rh/Re/Rn.
HAS_R_GROUP = re.compile(r"R(?![a-z])")

TOOL = "Marvin+RDKit"

STRUCT_ROOT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "Biochemistry", "Structures"
)

# Marvin prefixes its InChIKey export; the 23.4 bundles store the bare key.
INCHIKEY_PREFIX = re.compile(r"^InChIKey=")


def marvin_version(cxcalc="cxcalc"):
    """Major.minor of the installed suite, for the filename and tool_version column."""
    out = subprocess.run([cxcalc, "--help"], capture_output=True, text=True).stdout
    m = re.search(r"(\d+\.\d+)", out)
    return m.group(1) if m else "unknown"


def read_structures(source, column, filename):
    """Return [(external_id, structure)] in file order, skipping empty rows."""
    path = os.path.join(STRUCT_ROOT, source, filename)
    rows = []
    with open(path) as fh:
        for line in DictReader(fh, dialect="excel-tab"):
            ext_id, struct = line.get("external_id"), line.get(column)
            if ext_id and struct:
                rows.append((ext_id, struct))
    return rows


def _marvin_lib(cxcalc="cxcalc"):
    """The jar directory of the installed suite, for the plugin classpath."""
    lib = os.environ.get("MARVIN_LIB")
    if lib:
        return lib
    exe = shutil.which(cxcalc)
    if exe:
        cand = os.path.join(os.path.dirname(os.path.realpath(exe)), "..", "lib")
        if os.path.isdir(cand):
            return os.path.normpath(cand)
    raise RuntimeError(
        "cannot locate the Marvin jars; set MARVIN_LIB to <marvinsuite>/lib")


def _inchi_stereo_layer(inchi):
    """Stereo layers (/t, /m, /s, /b) of an InChI, or '' if none."""
    if not inchi or not inchi.startswith("InChI="):
        return ""
    parts = []
    for seg in inchi.split("/")[1:]:
        if seg and seg[0] in ("t", "m", "s", "b"):
            parts.append("/" + seg)
    return "".join(parts)


# Imidic acid / lactim O-H on a carbon double-bonded to nitrogen, matched on the
# Kekule form so that aromatic lactims (2-hydroxypyridine, the 4-hydroxy form of
# uracil, the 6-hydroxy form of guanine) are caught as well as acyclic ones.
_IMIDIC = None


def normalize_amide_lactam(rdmol, Chem):
    """Move H from an imidic-acid/lactim O to its N: -C(OH)=N-  ->  -C(=O)-NH-.

    WHY. RDKit's InChI reader draws many amides as imidic acids and many
    lactams as lactims, because InChI stores those H atoms as mobile. The
    canonical-tautomer pick repairs most of them, but not where the
    enumeration hits its cap (peptides, CoA esters) or where the stereo guard
    hands back the input. Measured on the bundles generated without this
    step: 686 compounds kept an amide N deprotonated as an imidate and 434 a
    lactam/guanine-type ring N-H, against amide N-H pKa ~15 and uridine N3
    9.25 / guanosine N1 9.25 measured (IUPAC). The same tautomers fed to the
    pKa run produce spurious pKas (uridine N3 6.28 from the InChI vs 9.15
    from the lactam).

    SAFETY. Applied only if the standard InChI of the result is IDENTICAL to
    the input's. Standard InChI treats these H as mobile, so an identical
    InChI proves the change is a pure tautomer move: no formula, charge,
    connectivity or stereo layer can change. Returns (mol, n_sites).
    """
    global _IMIDIC, _AMINO_OXO
    if _IMIDIC is None:
        _IMIDIC = Chem.MolFromSmarts("[O;X2;H1;+0]-[#6;X3;+0]=[#7;X2;+0]")
        # 2-amino-pyrimidin-4-one drawn with the ring H on the N away from the
        # C=O (RDKit's canonical pick for guanine and pterins). Guanine N1-H and
        # pterin N3-H -- the N flanked by C-NH2 and C=O -- is the measured form;
        # Marvin 26.1 gives guanosine 7.44 on RDKit's form vs 10.16 on N1-H
        # (measured 9.25).
        _AMINO_OXO = Chem.MolFromSmarts(
            "[NX3;!R]-c1:[nX2;H0;+0:1]:c(=O):c:c:[nX3;H1;+0:2]:1")
    try:
        inchi_in = Chem.MolToInchi(rdmol)
    except Exception:
        return rdmol, 0
    # Pass 1: 2-amino-pyrimidin-4-one ring H to the N between C-NH2 and C=O.
    n_ring = 0
    try:
        moved = Chem.RWMol(rdmol)
        for _ in range(20):
            hits = moved.GetSubstructMatches(_AMINO_OXO)
            if not hits:
                break
            m = hits[0]                        # [N, c2, n:1, c=, O, c, c, nH:2]
            to_n, from_n = moved.GetAtomWithIdx(m[2]), moved.GetAtomWithIdx(m[7])
            to_n.SetNumExplicitHs(1); to_n.SetNoImplicit(True)
            from_n.SetNumExplicitHs(0); from_n.SetNoImplicit(True)
            Chem.SanitizeMol(moved)
            n_ring += 1
        if n_ring:
            cand = moved.GetMol()
            if Chem.MolToInchi(cand) == inchi_in:
                rdmol = cand
            else:
                n_ring = 0
    except Exception:
        n_ring = 0

    # Pass 2: imidic acid / lactim O-H to N, on the Kekule form.
    n_imidic = 0
    try:
        work = Chem.RWMol(rdmol)
        Chem.Kekulize(work, clearAromaticFlags=True)
        for _ in range(200):                  # one site per pass; conjugation can expose more
            hits = work.GetSubstructMatches(_IMIDIC)
            if not hits:
                break
            o, c, nn = hits[0]
            work.GetBondBetweenAtoms(o, c).SetBondType(Chem.BondType.DOUBLE)
            work.GetBondBetweenAtoms(c, nn).SetBondType(Chem.BondType.SINGLE)
            ao, an = work.GetAtomWithIdx(o), work.GetAtomWithIdx(nn)
            h_n = an.GetTotalNumHs()
            ao.SetNumExplicitHs(0); ao.SetNoImplicit(True)
            an.SetNumExplicitHs(h_n + 1); an.SetNoImplicit(True)
            work.UpdatePropertyCache(strict=False)
            n_imidic += 1
        if n_imidic:
            out = work.GetMol()
            Chem.SanitizeMol(out)             # re-aromatizes lactam rings
            if Chem.MolToInchi(out) == inchi_in:
                return out, n_ring + n_imidic
    except Exception:
        pass
    return rdmol, n_ring


def make_tautomer_tools():
    """(TautomerEnumerator, Uncharger) configured as the shipped step uses them."""
    from rdkit.Chem.MolStandardize import rdMolStandardize as _rms
    te = _rms.TautomerEnumerator()
    te.SetMaxTautomers(50)
    te.SetMaxTransforms(1000)
    te.SetRemoveSp3Stereo(False)
    te.SetRemoveBondStereo(False)
    return te, _rms.Uncharger()


_PKA_TOOLS = None


def pka_input(struct):
    """Structure string Run_Marvin_pKas.py feeds Marvin: the same neutral,
    canonical, amide/lactam-normalized tautomer the protonation step uses,
    written as isomeric SMILES. None when RDKit cannot read the structure
    (the caller then feeds the original string). Override_Adenine_N1_pKas.py
    calls this too, so both see the identical atom order."""
    global _PKA_TOOLS
    from rdkit import Chem
    if _PKA_TOOLS is None:
        _PKA_TOOLS = make_tautomer_tools()
    try:
        mol = (Chem.MolFromInchi(struct) if struct.startswith("InChI=")
               else Chem.MolFromSmiles(struct))
        if mol is None:
            return None
        out, _ = tautomer_canonical(mol, _PKA_TOOLS[0], _PKA_TOOLS[1], Chem)
        smi = Chem.MolToSmiles(out)
        return smi or None
    except Exception:
        return None


def tautomer_canonical(rdmol, taut_enumerator, uncharger, Chem, notes=None):
    """Canonical tautomer, then amide/lactam normalization -- see _tautomer_core.

    notes, if a dict, receives notes["normalized"] = number of imidic-acid /
    lactim sites moved to amide / lactam by normalize_amide_lactam()."""
    mol, info = _tautomer_core(rdmol, taut_enumerator, uncharger, Chem)
    mol, n = normalize_amide_lactam(mol, Chem)
    if notes is not None:
        notes["normalized"] = n
    return mol, info


def _tautomer_core(rdmol, taut_enumerator, uncharger, Chem):
    """RDKit canonical tautomer of an RDKit Mol -- the (iii-RDKit) step.

    Replaces the Marvin Isomers Plugin Group TautomerizationPlugin
    (Step iii of ChargeMol.java), which requires a licence this
    regeneration does not hold.

    CHARGED SOURCES FIRST: MetaCyc/ChEBI/Rhea store charged species, and
    RDKit's InChI reader puts the negative charge on the oxygen, giving an
    imidate ANION N=C([O-]).  Tautomer enumeration never moves charge, so
    without neutralizing first the imidate defect survives untouched.  The
    input is therefore uncharged with rdMolStandardize.Uncharger before
    enumerating.  This is safe: Marvin's microspecies step recomputes the
    pH-7 charge afterwards, and the formula/charge gate still compares the
    final row against the source.

    CAP POLICY: the enumerator's MaxTautomers/MaxTransforms cap is the only
    wall-clock bound that works (a signal-based timeout cannot interrupt
    the C++ enumeration).  When the cap is hit the (partial) canonical pick
    is STILL USED and the cap counted as "tautomer_cap" -- abandoning the
    compound leaves the defect in place; peptides and CoA esters are
    exactly the molecules that hit the cap.

    STEREO GUARD: the InChI stereo layers of input vs output are compared;
    a loss of any /t,/m,/s,/b layer returns the (uncharged) input instead
    -- stereochemistry takes precedence over tautomer normalization, and the
    downstream gate() only checks formula/charge and would not catch it.

    Returns (mol, info) where info is one of "ok", "tautomer_cap",
    "stereo_loss", "fallback".
    """
    if taut_enumerator is None:
        return rdmol, "fallback"

    try:
        inchi_before = Chem.MolToInchi(rdmol) if rdmol else ""
    except Exception:
        inchi_before = ""
    stereo_before = _inchi_stereo_layer(inchi_before)

    rdmol_n = uncharger.uncharge(rdmol) if uncharger is not None else rdmol

    try:
        res = taut_enumerator.Enumerate(rdmol_n)
        capped = False
        try:
            if getattr(res, "status", None) is not None:
                from rdkit.Chem.MolStandardize.rdMolStandardize import (
                    TautomerEnumeratorStatus)
                capped = res.status != TautomerEnumeratorStatus.Completed
        except Exception:
            pass
        canon = taut_enumerator.PickCanonical(res)
        if canon is None:
            return rdmol_n, "fallback"
    except Exception:
        return rdmol_n, "fallback"

    if stereo_before:
        try:
            stereo_after = _inchi_stereo_layer(Chem.MolToInchi(canon) if canon else "")
            if stereo_before != stereo_after:
                return rdmol_n, "stereo_loss"
        except Exception:
            pass

    return canon, ("tautomer_cap" if capped else "ok")


class Protonator:
    """MajorMicrospeciesPlugin at a fixed pH, one fragment at a time.

    Mirrors ChargeMol.java minus its tautomer step -- see METHODOLOGY DELTA.
    The JVM starts on construction, so --help never pays for it.
    """

    INCHI_FORMAT = "inchi:AuxNone,SAbs"        # SAbs is load-bearing -- see EXPORT

    def __init__(self, ph=7.0, cxcalc="cxcalc"):
        import jpype                                   # noqa: local by design
        from rdkit import Chem, RDLogger               # noqa: the import bridge
        RDLogger.DisableLog("rdApp.*")                 # parse notes are not findings
        self._Chem = Chem
        # RDKit tautomer enumerator — replaces Marvin Isomers TautomerizationPlugin
        # (Step iii of ChargeMol.java).  See _tautomer_canonical.
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize as _rms
            self._taut_enumerator, self._uncharger = make_tautomer_tools()
            self._rdMolStandardize = _rms
        except Exception:
            self._taut_enumerator = None
            self._uncharger = None
            self._rdMolStandardize = None
        if not jpype.isJVMStarted():
            jpype.startJVM(classpath=glob.glob(os.path.join(_marvin_lib(cxcalc), "*.jar")))
        self._jpype = jpype
        self.Molecule = jpype.JClass("chemaxon.struc.Molecule")
        self.MolImporter = jpype.JClass("chemaxon.formats.MolImporter")
        self.MolExporter = jpype.JClass("chemaxon.formats.MolExporter")
        self.plugin = jpype.JClass("chemaxon.calculations.MajorMicrospeciesPlugin")()
        self.plugin.setpH(float(ph))

    @staticmethod
    def _inchi_stereo_layer(inchi):
        return _inchi_stereo_layer(inchi)

    def _tautomer_canonical(self, rdmol, notes=None):
        """See the module-level tautomer_canonical() for the full contract.

        Uncharges first (charged source InChIs parse as imidate anions,
        which tautomer enumeration cannot repair), uses the canonical
        pick even when the enumerator cap is hit (counted as
        "tautomer_cap", not abandoned), and falls back to the uncharged
        input on stereo loss.
        """
        return tautomer_canonical(
            rdmol, self._taut_enumerator, self._uncharger, self._Chem, notes=notes)

    # Per-compound tautomer state: dedup counter inflation from
    # multiple calls (aromatized/unaromatized × InChI/SMILES ×
    # disconnected-InChI row).  process() clears this each compound.
    _seen_tautomer_ids = None

    def _import(self, struct, stats=None):
        """Marvin Molecule for a structure string, bridged through RDKit.

        Marvin's own InChI and SMILES importers drop stereochemistry (61% and
        39% round-trip fidelity); its molblock importer does not (99.6%). So
        RDKit parses the string, uncharges it (charged sources parse as
        imidate anions, which tautomer enumeration cannot repair), picks
        its canonical tautomer (replacing Marvin's Isomers-licensed step
        iii), and hands Marvin a molblock instead. Falls back to a direct
        import only when RDKit cannot read the structure at all, so that an
        odd structure still gets a chance rather than being dropped.
        stats, if given, collects tautomer counters (cap hits, stereo
        loss, fallbacks) once per compound.
        """
        Chem = self._Chem
        rdmol = (Chem.MolFromInchi(struct) if struct.startswith("InChI=")
                 else Chem.MolFromSmiles(struct))
        if rdmol is not None:
            notes = {}
            rdmol, info = self._tautomer_canonical(rdmol, notes=notes)
            if stats is not None and notes.get("normalized"):
                seen = self._seen_tautomer_ids
                if seen is None:
                    seen = self._seen_tautomer_ids = set()
                key = (stats.get("_current_id"), "amide_lactam")
                if key not in seen:
                    seen.add(key)
                    stats["amide_lactam_normalized"] = stats.get("amide_lactam_normalized", 0) + 1
            # stats dedup: one (external_id, struct) may be tried up to
            # four times (aromatized/unaromatized × candidates); count
            # each external_id once per non-ok outcome.
            if stats is not None and info != "ok":
                seen = self._seen_tautomer_ids
                if seen is None:
                    seen = self._seen_tautomer_ids = set()
                key = (stats.get("_current_id"), info)
                if key not in seen:
                    seen.add(key)
                    if info == "tautomer_cap":
                        stats["tautomer_cap_hit"] = stats.get("tautomer_cap_hit", 0) + 1
                    elif info == "stereo_loss":
                        stats["tautomer_stereo_loss"] = stats.get("tautomer_stereo_loss", 0) + 1
                    elif info == "fallback":
                        stats["tautomer_fallback"] = stats.get("tautomer_fallback", 0) + 1
            block = Chem.MolToMolBlock(rdmol, kekulize=True)
            if block:
                return self.MolImporter.importMol(block)
        return self.MolImporter.importMol(struct)

    def fragment_count(self, struct):
        """Disconnected components RDKit sees in a structure string, or None.

        Used to catch metal disconnection -- see INPUT SELECTION.
        """
        Chem = self._Chem
        try:
            mol = (Chem.MolFromInchi(struct) if struct.startswith("InChI=")
                   else Chem.MolFromSmiles(struct))
            return len(Chem.GetMolFrags(mol)) if mol is not None else None
        except Exception:
            return None

    def protonate(self, struct, aromatize=True, stats=None):
        """Protonated Molecule for one structure string, or None if unusable."""
        mol = self._import(struct, stats=stats)
        fused = self.Molecule()
        for frag in mol.convertToFrags():
            if aromatize:
                frag.aromatize()
            self.plugin.setMolecule(frag)
            self.plugin.run()
            fused.fuse(self.plugin.getMajorMicrospecies(), False)
        if fused.getAtomCount() == 0:
            return None
        # Adenine N1: 26.1 puts its pKb at ~7.0 and protonates it; 23.4 put
        # it at ~2.4. Use 23.4's -- see Adenine_N1.py and ADENINE N1 above.
        if Adenine_N1.deprotonate_n1(fused, self.MolExporter) and stats is not None:
            seen = self._seen_tautomer_ids
            if seen is None:
                seen = self._seen_tautomer_ids = set()
            key = (stats.get("_current_id"), "adenine_n1")
            if key not in seen:
                seen.add(key)
                stats["adenine_n1_override"] = stats.get("adenine_n1_override", 0) + 1
        return fused

    def protonate_best(self, candidates, stats=None):
        """First (molecule, smiles) from a ladder of (structure, aromatize) tries.

        A protonated molecule is only useful if it can be written back out, and
        two things make that fail on a small tail of compounds:

          aromatize() on a QUERY molecule leaves bonds no SMILES writer can
          express, so the export throws and the compound would be dropped. On
          ChEBI 26432 the un-aromatized retry writes `*C1=[NH+]C(=*)NC(*)=C1*`
          -- character for character what 23.4 shipped.

          Some InChIs do not round-trip to a writable molecule even though the
          compound's own SMILES does. Azide is the clean example: the source
          SMILES `[N-]=[N+]=[N-]` protonates and exports fine, while its InChI
          (which does not carry the charge distribution) rebuilds as a radical
          that the SMILES writer rejects.

        So each candidate is tried aromatized first -- ChargeMol.java's
        behaviour, and what the bulk of the file is generated with -- then
        un-aromatized, before falling back to the next representation. The
        ladder only ever advances when the current rung produces nothing
        writable, so it cannot change a compound that already worked.
        """
        rungs = [(struct, aromatize)
                 for struct in candidates for aromatize in (True, False)]

        # Pass 1 -- Marvin's own SMILES writer on every rung. Its refusal is the
        # signal that drives the ladder, so the whole ladder has to be walked
        # with it before any rescue is attempted.
        molecules = []
        for struct, aromatize in rungs:
            try:
                mol = self.protonate(struct, aromatize=aromatize, stats=stats)
            except Exception:
                continue
            if mol is None:
                continue
            molecules.append(mol)
            smiles = self.export(mol, "smiles")
            if smiles:
                return mol, smiles

        # Pass 2 -- only now, with every rung refused, fall back to the molblock
        # route for the query-bond structures Marvin cannot write at all. Doing
        # this inside pass 1 would let the first rung's degraded output win over
        # a later rung's clean one: on ChEBI 26432 the aromatized rung rescues
        # to `*C1~N~C(=*)~[N+]~C(*)~C~1*` while the un-aromatized rung writes
        # 23.4's `*C1=[NH+]C(=*)NC(*)=C1*`.
        for mol in molecules:
            smiles = self.smiles_of(mol)
            if smiles:
                return mol, smiles
        return None, ""

    def export(self, mol, fmt):
        """First line of an export, or '' -- Marvin appends AuxInfo to some formats."""
        try:
            out = str(self.MolExporter.exportToFormat(mol, fmt)).strip()
        except Exception:
            return ""
        return out.split("\n")[0].strip() if out else ""

    def inchikey(self, inchi):
        """The InChIKey OF a given InChI string -- see EXPORT."""
        try:
            return INCHIKEY_PREFIX.sub("", self._Chem.InchiToInchiKey(inchi) or "")
        except Exception:
            return ""

    def smiles_of(self, mol):
        """SMILES for a protonated molecule, Marvin's writer first.

        Marvin cannot write plain SMILES for a molecule carrying QUERY bonds --
        MetaCyc spells its protein-bound cofactors with `~` (any-bond), as in
        `*~O=C/C=C(C)/...` for the retinal-binding proteins, and the writer
        throws rather than degrading them. Going out through a molblock and
        letting RDKit write the SMILES keeps the query bond: on
        11C-Retinal-RALBPs that route returns the 23.4 string character for
        character. sanitize=False because these structures are deliberately
        not valid closed-shell molecules and sanitizing would reject them.
        """
        direct = self.export(mol, "smiles")
        if direct:
            return direct
        try:
            block = str(self.MolExporter.exportToFormat(mol, "mol"))
            rdmol = self._Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)
            return self._Chem.MolToSmiles(rdmol) if rdmol is not None else ""
        except Exception:
            return ""


def process(source, version, ph, limit=0, cxcalc="cxcalc"):
    """Protonate one source and write its bundle. Returns (path, stats)."""
    inchi_rows = read_structures(source, "inchi", "inchi.tsv")
    smiles_rows = read_structures(source, "smiles", "smiles.tsv")
    inchi_by_id = dict(inchi_rows)

    # smiles.tsv is the superset, so it defines the compound list and the row
    # order. InChI wins as the INPUT wherever it exists -- see INPUT SELECTION.
    compounds = smiles_rows[:limit] if limit else smiles_rows

    prot = Protonator(ph=ph, cxcalc=cxcalc)
    generated_on = date.today().isoformat()
    # The 23.4 bundles store a bare "7"; don't drift to "7.0".
    ph_out = int(ph) if float(ph).is_integer() else ph

    out_rows = []
    # Failure taxonomy, kept separate so "no result" is never reported as though
    # it were all failure -- the pKa run had to walk that conflation back once.
    stats = {"in": len(compounds), "ok": 0, "unparseable": 0,
             "plugin_error": 0, "empty": 0, "nonstandard_inchi": 0,
             "formula_from_marvin": 0, "wildcard_r_added": 0,
             "smiles_preferred_disconnected_inchi": 0, "inchi_from_inchi_form": 0,
             "passthrough_invariant": 0,
             "tautomer_fallback": 0, "tautomer_cap_hit": 0,
             "tautomer_stereo_loss": 0, "adenine_n1_override": 0,
             "amide_lactam_normalized": 0}
    # Source formula/charge per (external_id, type), the reference every
    # protonated row is gated against. inchi.tsv/smiles.tsv, already on disk.
    _inv.STRUCT = STRUCT_ROOT
    source_rows = _inv.load_source(source)
    passed_through = []

    def gate(struct_type, structure, formula, charge):
        """Return (structure, formula, charge) to write for this row: the
        protonated triple when it obeys the invariant, the SOURCE triple when
        it does not. INFO kinds (no source formula to compare against, or a
        wildcard structure whose H count is convention-dependent) are never
        replaced -- the rule only touches rows the invariant rejects."""
        src = source_rows.get((ext_id, struct_type))
        if not src:
            return structure, formula, charge
        kind, dH, dq = _inv.compare(src[0], src[1], formula, charge)
        if kind == "ok" or kind in _inv.INFO_KINDS:
            return structure, formula, charge
        stats["passthrough_invariant"] += 1
        passed_through.append((ext_id, struct_type, kind, dH, dq,
                               structure, formula, charge, src[2], src[0], src[1]))
        return src[2], src[0], src[1]

    for ext_id, smiles in compounds:
        stats["_current_id"] = ext_id
        prot._seen_tautomer_ids = set()
        # InChI first, the compound's own SMILES as the fallback rung -- unless
        # the InChI is the MORE fragmented of the two, which means it has
        # disconnected something the SMILES keeps bonded. See INPUT SELECTION.
        last_resort = []
        disconnected = False
        if ext_id in inchi_by_id:
            inchi = inchi_by_id[ext_id]
            f_inchi = prot.fragment_count(inchi)
            f_smiles = prot.fragment_count(smiles)
            # RDKit failing to read the SMILES is not evidence it is the worse
            # input -- Marvin reads every one of these clusters fine. Treat an
            # unknown SMILES count as "not more fragmented".
            disconnected = (f_inchi is not None and f_inchi > 1
                            and (f_smiles is None or f_smiles < f_inchi))
            if disconnected:
                # DEMOTE the InChI rather than just reordering. Sharing one
                # ladder lets a writable-but-wrong rung beat a correct one: on
                # KEGG C18384 the SMILES yields the right dative-bonded
                # magnesium propionate, which Marvin's SMILES writer refuses,
                # while the disconnected InChI yields `[Mg++].CCC([O-])=O.
                # CCC([O-])=O` -- writable, and three fragments where 23.4 has
                # one. The InChI stays only as a last resort, so it is reached
                # when the SMILES gives NOTHING rather than when it merely
                # needs the molblock rescue.
                stats["smiles_preferred_disconnected_inchi"] += 1
                candidates, last_resort = [smiles], [inchi]
            else:
                candidates = [inchi, smiles]
        else:
            candidates = [smiles]
        try:
            mol, smile_out = prot.protonate_best(candidates, stats=stats)
            if (mol is None or not smile_out) and last_resort:
                mol, smile_out = prot.protonate_best(last_resort, stats=stats)
        except Exception as exc:
            name = type(exc).__name__
            # A rejected structure and a plugin crash are different findings.
            if "MolFormat" in name or "Import" in name:
                stats["unparseable"] += 1
            else:
                stats["plugin_error"] += 1
            continue
        if mol is None or not smile_out:
            # Nothing writable came out of any rung. For H+ (InChI=1S/p+1) that
            # is the right answer and 23.4 agrees -- it has no row either.
            stats["empty"] += 1
            continue

        def columns(struct_type, structure):
            """(formula, charge) for one row, the way this repository derives them.

            Per ROW, not per compound, because refresh_file() recomputes each
            row from its own structure string and the two representations can
            disagree. Marvin's own getFormula() is the fallback only.
            """
            try:
                f, c, _ = parse_structure(struct_type, structure)
            except Exception:
                f = c = None
            if f is None:
                # Neither RDKit nor OpenBabel could read it -- typically a
                # deliberately invalid valence like ISOCITHASE-P's
                # `*OP(=O)(=O)=O`. Marvin's formula is the only one available.
                stats["formula_from_marvin"] += 1
                f, c = str(mol.getFormula()), str(mol.getTotalCharge())

            # INVARIANT: a structure carrying `*` gets an R in its formula.
            # parse_structure implements that as re.sub(r'\\*', 'R', formula),
            # which only works when RDKit produced the formula -- RDKit renders
            # a dummy atom as `*`, OpenBabel and Marvin both omit it entirely.
            # So on the OpenBabel path the substitution finds nothing to rewrite
            # and the row ships a formula contradicting its own structure
            # column. Enforcing the invariant here completes the convention; it
            # does not re-derive the formula, and it is asserted over the whole
            # bundle afterwards so it cannot regress silently again.
            if not HAS_R_GROUP.search(f):
                wildcards = sum(
                    1 for i in range(mol.getAtomCount())
                    if str(mol.getAtom(i).getSymbol()) in WILDCARD_SYMBOLS)
                if wildcards:
                    stats["wildcard_r_added"] += 1
                    f = Compounds.mergeFormula(
                        f + ("R" if wildcards == 1 else f"R{wildcards}"))[0]
            return f, c

        stats["ok"] += 1
        smile_formula, smile_charge = columns("SMILE", smile_out)
        smile_out, smile_formula, smile_charge = gate(
            "SMILE", smile_out, smile_formula, smile_charge)
        ver_str = _rdk_ver_str(version)
        out_rows.append((ext_id, "SMILE", smile_out, smile_formula, smile_charge,
                         TOOL, ver_str, ph_out, generated_on))

        # InChI and InChIKey only where the source carries an InChI: InChI
        # cannot represent the query molecules -- see OUTPUT LAYOUT.
        if ext_id in inchi_by_id:
            # The InChI row is written from the INCHI-derived molecule when the
            # two representations disagree about connectivity, because each
            # column should carry what its own representation can express --
            # which is exactly what 23.4 did. For CPD-18407 it shipped a
            # connected 8-iron cluster in SMILE and the disconnected
            # `InChI=1S/C.8Fe.6HS.3S/...` in InChI.
            #
            # This is also load-bearing for stability, not just fidelity.
            # Asking Marvin to write an InChI for a CONNECTED metal cluster
            # aborts the JVM outright: `free(): double free detected in tcache
            # 2` from InChINativeGenerateInChICall, a native fault no Python or
            # Java handler can catch. Feeding the InChI writer the structure
            # InChI can actually represent avoids the crash by construction.
            # Protonate the INCHI-derived molecule for this row when the two
            # representations disagree about connectivity. Each column then
            # carries what its own representation can express, which is what
            # 23.4 did -- CPD-18407 ships a connected 8-iron cluster in SMILE
            # and the disconnected `InChI=1S/C.8Fe.6HS.3S/...` in InChI.
            #
            # Protonate it rather than passing the source through. 23.4
            # protonated this form uniformly; the protonation merely happened
            # to be a no-op for 234 of the 493, and was a real change for the
            # other 259 (chlorophylls and cobalamins pick up a /p-2 layer).
            # Passing the source through would silently un-protonate those 259.
            #
            # This is also load bearing for stability. Asking Marvin to write
            # an InChI for a CONNECTED metal cluster aborts the JVM outright:
            # `free(): double free detected in tcache 2` from
            # InChINativeGenerateInChICall, a native fault no Python or Java
            # handler can catch, and one a previous regeneration died on.
            # Feeding the InChI writer the structure InChI can actually
            # represent avoids the crash by construction.
            inchi_source = mol
            if disconnected:
                stats["inchi_from_inchi_form"] += 1
                inchi_source, _ = prot.protonate_best([inchi_by_id[ext_id]], stats=stats)
            inchi_out = (prot.export(inchi_source, Protonator.INCHI_FORMAT)
                         if inchi_source is not None else "")
            if inchi_out:
                if not inchi_out.startswith("InChI=1S/"):
                    # Standard InChI only; a non-standard string would not be
                    # comparable with the rest of the column.
                    stats["nonstandard_inchi"] += 1
                else:
                    inchi_formula, inchi_charge = columns("InChI", inchi_out)
                    inchi_out, inchi_formula, inchi_charge = gate(
                        "InChI", inchi_out, inchi_formula, inchi_charge)
                    out_rows.append((ext_id, "InChI", inchi_out,
                                     inchi_formula, inchi_charge,
                                     TOOL, ver_str, ph_out, generated_on))
                    # Hash the string just written, never a second export -- see
                    # EXPORT for why Marvin's own inchikey disagrees with it.
                    # After the gate that string may be the SOURCE InChI, and
                    # the key must follow it.
                    key_out = prot.inchikey(inchi_out)
                    if key_out:
                        out_rows.append((ext_id, "InChIKey", key_out, "", "",
                                         TOOL, ver_str, ph_out, generated_on))

    ph_tag = f"ph{int(ph)}" if float(ph).is_integer() else f"ph{ph}"
    out_dir = os.path.join(STRUCT_ROOT, source, "protonations")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"marvin_{version}_{ph_tag}.tsv")
    with open(out_path, "w") as fh:
        fh.write("external_id\ttype\tstructure\tformula\tcharge\t"
                 "tool\ttool_version\tph\tgenerated_on\n")
        for row in out_rows:
            fh.write("\t".join(str(c) for c in row) + "\n")

    # Every row the gate replaced, with what it replaced -- next to the other
    # run reports, NOT under protonations/, which BiochemPy globs as bundles.
    if passed_through:
        rep_dir = os.path.join(STRUCT_ROOT, "_reports")
        os.makedirs(rep_dir, exist_ok=True)
        side = os.path.join(rep_dir, f"marvin_{version}_{ph_tag}_passthrough_{source}.tsv")
        with open(side, "w") as fh:
            fh.write("external_id\ttype\tkind\tdH\tdq\twas_structure\twas_formula\t"
                     "was_charge\tnow_structure\tnow_formula\tnow_charge\n")
            for row in passed_through:
                fh.write("\t".join(str(c) for c in row) + "\n")

    # _current_id is internal loop state, not a reported stat
    stats.pop("_current_id", None)
    print(f"{source:<8} in={stats['in']:<6} protonated={stats['ok']:<6} "
          f"unparseable={stats['unparseable']:<4} plugin_error={stats['plugin_error']:<4} "
          f"empty={stats['empty']:<4} nonstandard_inchi={stats['nonstandard_inchi']:<4} "
          f"formula_from_marvin={stats['formula_from_marvin']:<4} "
          f"wildcard_r_added={stats['wildcard_r_added']:<4} "
          f"smiles_preferred={stats['smiles_preferred_disconnected_inchi']:<4} "
          f"inchi_from_inchi_form={stats['inchi_from_inchi_form']:<4} "
          f"passthrough_invariant={stats['passthrough_invariant']} "
          f"tautomer_fallback={stats['tautomer_fallback']:<4} "
          f"tautomer_cap={stats['tautomer_cap_hit']:<4} "
          f"tautomer_stereo_loss={stats['tautomer_stereo_loss']:<4} "
          f"adenine_n1_override={stats['adenine_n1_override']:<4} "
          f"amide_lactam_normalized={stats['amide_lactam_normalized']:<4}")
    print(f"{'':<8} rows={len(out_rows):<7} -> {os.path.relpath(out_path, STRUCT_ROOT)}")
    return out_path, stats


def main():
    cxcalc = os.environ.get("CXCALC_BIN", "cxcalc")
    version = marvin_version(cxcalc)
    try:
        from rdkit import __version__ as _rdk_ver
    except Exception:
        _rdk_ver = "unknown"
    print(f"Marvin {version} MajorMicrospeciesPlugin + RDKit {_rdk_ver} "
          f"TautomerEnumerator (tautomer step iii-RDKit) + adenine N1 at 23.4's pKa, "
          f"pH {_ARGS.ph}\n")
    for source in (_ARGS.sources or SOURCES):
        process(source, version, _ARGS.ph, limit=_ARGS.limit, cxcalc=cxcalc)


if __name__ == "__main__":
    main()
