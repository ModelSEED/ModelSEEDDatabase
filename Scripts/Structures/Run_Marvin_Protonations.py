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
    iii) find the dominant tautomer at pH 7   <-- NOT REPRODUCED HERE
    iv)  find the major microspecies at pH 7
    v)   convert implicit hydrogens back to explicit

Step (iii) needs chemaxon TautomerizationPlugin, which is licensed under the
ISOMERS Plugin Group. The license this repository's regeneration runs on covers
the PROTONATION Plugin Group only -- pKa and major microspecies. The refusal is
explicit and was verified, not assumed:

    $ cxcalc "CC(=O)CC(=O)C" majortautomer -H 7
    chemaxon.license.api.LicenseException: No valid license has been found.
    Product name: Isomers Plugin Group

So this script runs (i), (ii), (iv), (v) and SKIPS (iii). The consequence is
measurable and is not small. Across the 53,021 compounds both bundles cover,
the net charge at pH 7 is unchanged for 79.7% and different for 20.3%, skewed
toward deprotonation: 11.8% of compounds sit one charge unit lower under 26.1
and 4.3% two or more lower, against 3.6% one higher and 0.7% two or more.

That 20.3% is an UPPER BOUND on the tautomer effect, not a measurement of it.
It also contains the genuine 23.4 -> 26.1 engine improvement, which is already
known to be substantial: the pKa regeneration found 26.1 disagreeing with 23.4
about the NUMBER of ionizable sites in roughly one shared set in seven. A pKa
that moves across 7 flips the protonation state at pH 7 on its own, with no
tautomer involved. The two causes cannot be separated without an Isomers
license, and this script does not pretend otherwise.

If that license is obtained, restore step (iii) by inserting a
TautomerizationPlugin call (setpH(7.0), setDominantTautomerDistribution-
Calculation(true), dominant tautomer is getStructure(0)) between the aromatize
and the microspecies call in protonate(), and regenerate. Nothing else changes.

INPUT SELECTION follows Run_Marvin_pKas.py: `smiles.tsv` is a strict superset
of `inchi.tsv`, and InChI wins wherever it exists, because a charged SMILES is
an already-deprotonated species and Marvin would be answering a different
question than it does for the neutral InChI parent. Each compound is protonated
exactly once, from its preferred representation.

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

# What Marvin calls a wildcard atom, across the notations it accepts. `R#` is
# the one that matters and the one easiest to miss: a molblock R atom -- which
# is how RDKit writes every dummy atom, and therefore how EVERY structure
# arrives here through the import bridge -- reads back as symbol "R#", not "R".
# A set without it counts zero wildcards on a molecule that plainly has them.
WILDCARD_SYMBOLS = {"*", "A", "R", "R#"}

# An R group in a formula, not the R of Ru/Rb/Rh/Re/Rn.
HAS_R_GROUP = re.compile(r"R(?![a-z])")

TOOL = "Marvin"

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
        if not jpype.isJVMStarted():
            jpype.startJVM(classpath=glob.glob(os.path.join(_marvin_lib(cxcalc), "*.jar")))
        self._jpype = jpype
        self.Molecule = jpype.JClass("chemaxon.struc.Molecule")
        self.MolImporter = jpype.JClass("chemaxon.formats.MolImporter")
        self.MolExporter = jpype.JClass("chemaxon.formats.MolExporter")
        self.plugin = jpype.JClass("chemaxon.calculations.MajorMicrospeciesPlugin")()
        self.plugin.setpH(float(ph))

    def _import(self, struct):
        """Marvin Molecule for a structure string, bridged through RDKit.

        Marvin's own InChI and SMILES importers drop stereochemistry (61% and
        39% round-trip fidelity); its molblock importer does not (99.6%). So
        RDKit parses the string and hands Marvin a molblock instead. Falls back
        to a direct import only when RDKit cannot read the structure at all, so
        that an odd structure still gets a chance rather than being dropped.
        """
        Chem = self._Chem
        rdmol = (Chem.MolFromInchi(struct) if struct.startswith("InChI=")
                 else Chem.MolFromSmiles(struct))
        if rdmol is not None:
            block = Chem.MolToMolBlock(rdmol, kekulize=True)
            if block:
                return self.MolImporter.importMol(block)
        return self.MolImporter.importMol(struct)

    def protonate(self, struct, aromatize=True):
        """Protonated Molecule for one structure string, or None if unusable."""
        mol = self._import(struct)
        fused = self.Molecule()
        for frag in mol.convertToFrags():
            if aromatize:
                frag.aromatize()
            self.plugin.setMolecule(frag)
            self.plugin.run()
            fused.fuse(self.plugin.getMajorMicrospecies(), False)
        if fused.getAtomCount() == 0:
            return None
        return fused

    def protonate_best(self, candidates):
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
                mol = self.protonate(struct, aromatize=aromatize)
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
             "formula_from_marvin": 0, "wildcard_r_added": 0}

    for ext_id, smiles in compounds:
        # InChI first, the compound's own SMILES as the fallback rung.
        candidates = [inchi_by_id[ext_id], smiles] if ext_id in inchi_by_id else [smiles]
        try:
            mol, smile_out = prot.protonate_best(candidates)
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
        out_rows.append((ext_id, "SMILE", smile_out, smile_formula, smile_charge,
                         TOOL, version, ph_out, generated_on))

        # InChI and InChIKey only where the source carries an InChI: InChI
        # cannot represent the query molecules -- see OUTPUT LAYOUT.
        if ext_id in inchi_by_id:
            inchi_out = prot.export(mol, Protonator.INCHI_FORMAT)
            if inchi_out:
                if not inchi_out.startswith("InChI=1S/"):
                    # Standard InChI only; a non-standard string would not be
                    # comparable with the rest of the column.
                    stats["nonstandard_inchi"] += 1
                else:
                    inchi_formula, inchi_charge = columns("InChI", inchi_out)
                    out_rows.append((ext_id, "InChI", inchi_out,
                                     inchi_formula, inchi_charge,
                                     TOOL, version, ph_out, generated_on))
                    # Hash the string just written, never a second export -- see
                    # EXPORT for why Marvin's own inchikey disagrees with it.
                    key_out = prot.inchikey(inchi_out)
                    if key_out:
                        out_rows.append((ext_id, "InChIKey", key_out, "", "",
                                         TOOL, version, ph_out, generated_on))

    ph_tag = f"ph{int(ph)}" if float(ph).is_integer() else f"ph{ph}"
    out_dir = os.path.join(STRUCT_ROOT, source, "protonations")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"marvin_{version}_{ph_tag}.tsv")
    with open(out_path, "w") as fh:
        fh.write("external_id\ttype\tstructure\tformula\tcharge\t"
                 "tool\ttool_version\tph\tgenerated_on\n")
        for row in out_rows:
            fh.write("\t".join(str(c) for c in row) + "\n")

    print(f"{source:<8} in={stats['in']:<6} protonated={stats['ok']:<6} "
          f"unparseable={stats['unparseable']:<4} plugin_error={stats['plugin_error']:<4} "
          f"empty={stats['empty']:<4} nonstandard_inchi={stats['nonstandard_inchi']:<4} "
          f"formula_from_marvin={stats['formula_from_marvin']:<4} "
          f"wildcard_r_added={stats['wildcard_r_added']}")
    print(f"{'':<8} rows={len(out_rows):<7} -> {os.path.relpath(out_path, STRUCT_ROOT)}")
    return out_path, stats


def main():
    cxcalc = os.environ.get("CXCALC_BIN", "cxcalc")
    version = marvin_version(cxcalc)
    print(f"Marvin {version} MajorMicrospeciesPlugin, pH {_ARGS.ph} "
          f"(no tautomer step -- see module docstring)\n")
    for source in (_ARGS.sources or SOURCES):
        process(source, version, _ARGS.ph, limit=_ARGS.limit, cxcalc=cxcalc)


if __name__ == "__main__":
    main()
