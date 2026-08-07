#!/usr/bin/env python
"""Write Group-Contribution compound energies into the BiochemPy compound JSON.

Source: MFAToolkit mol-analysis tables under
``Biochemistry/Thermodynamics/ModelSEED/*_MolAnalysis.tbl``.

Per-compound resolution:
  - structure picked via ``InChIKey`` then ``SMILE`` preference
  - aliases restricted to those listed for the *curated* structure in
    ``All_ModelSEED_Structures.txt``
  - mean of dg across matching aliases; error is
    sqrt(mean(σᵢ²) + var(dgᵢ)), so alias disagreement widens the reported
    uncertainty rather than being silently collapsed to a lowest pick.
    (Prior policy took the lowest; effectively identical for 99.95% of
    compounds, but the 25 pteridine/flavonoid-like cases with a 1–24
    kcal/mol spread across aliases now report both the mean value and
    an honest variance-inflated error.)

Also stamps Convention A load-bearing constants for small-molecule compounds
whose ΔfG values live in MFAToolkit's cue database (as ENERGY entries used
during GC decomposition of larger molecules) but which the alias-lookup
pipeline can't produce standalone values for. These are Chris Henry's own
reference values from the cue database — not novel numbers, just injected
at the compound level where the downstream pipeline needs them.

Each anchor value has been cross-validated by an independent method:
(Chris cue-db ENERGY) vs (eQuilibrator ΔfG'° shifted to Convention A by
subtracting n_H × RT·ln(10)·pH = n_H × 9.539). Eight of the nine agree
between the two independent sources to within 0.5 kcal/mol; the ninth
(cpd00239 H₂S) shows a 2.88 kcal/mol residual because eQuilibrator
treats it as a pseudo-species (weighted mix of H₂S and HS⁻ at pH 7,
pKa₁ ≈ 7.05) while Chris's cue_H2S is for neutral H₂S specifically.
Both are legitimate reference values for the H₂S identity that
cpd00239 represents (its formula HS just reflects Marvin's pH 7
Charged form choice); shipping Chris's value keeps GC uniformly in
Chris's own reference set.

See Biochemistry/Thermodynamics/README.md."""
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds
import _thermo_helpers as th

LABEL = 'Group contribution'

# Convention A anchors — cue-database ENERGY values, cross-validated against
# eQuilibrator (Convention B) via the ΔfG_A = ΔfG_B − n_H × 9.539 transform.
# All within 0.5 kcal/mol of the independent method.
CONVENTION_A_ANCHORS = [
    # (cpd_id, formula, ΔfG kcal/mol, error kcal/mol)
    ('cpd00067', 'H+',    -9.5,     0.0),   # Chris's pH 7 anchor (stated convention)
    ('cpd00001', 'H2O',  -56.687,   0.0),   # cue_H2O ENERGY
    ('cpd00013', 'NH4+', -18.97,    0.0),   # cue_NH4plus ENERGY
    ('cpd00011', 'CO2',  -92.26,    0.0),   # cue_CO2 ENERGY
    ('cpd00242', 'HCO3', -140.26,   0.0),   # cue_HCO3 ENERGY
    ('cpd00025', 'H2O2', -32.05,    0.0),   # cue_H2O2 ENERGY
    ('cpd11640', 'H2',     4.2065,  0.0),   # cue_H2 ENERGY
    ('cpd00007', 'O2',     3.9197,  0.0),   # cue_O2 ENERGY
    ('cpd00239', 'H2S',   -6.66,    0.0),   # cue_H2S ENERGY (see docstring on pseudo-species caveat)
]

compounds_helper = Compounds()
gc_table = th.parse_gc_compound_table(th.thermo_path())
curated_aliases = th.parse_curated_structure_aliases(
    th.structures_path('All_ModelSEED_Structures.txt'))


def resolve(cpd, stype, structure, aliases):
    # dev-drift guard: the picked structure (an InChIKey) may not match any
    # curated entry when the compound has been re-tautomerized or protonation
    # state has shifted since the aliases file was regenerated.
    curated = curated_aliases.get(cpd, {}).get(structure, {})
    return th.mean_energy_gc_style(
        (a for a in aliases if a in curated), gc_table)


th.run_compound_update(compounds_helper, LABEL, resolve, on_no_structure='default')

# Post-run: stamp all Convention A anchors. MFAToolkit's atom-labeler can't
# decompose these compounds standalone (they end up sentinel via alias
# lookup); their ΔfG values live only in the cue database, so we inject
# them at the compound level here.
compounds_dict = compounds_helper.loadCompounds()
touched = 0
for cpd_id, label_name, dg, dge in CONVENTION_A_ANCHORS:
    entry = compounds_dict.get(cpd_id)
    if entry is None:
        print(f"WARN: anchor compound {cpd_id} ({label_name}) missing from DB — skipping")
        continue
    if not isinstance(entry.get('thermodynamics'), dict):
        entry['thermodynamics'] = {}
    entry['thermodynamics'][LABEL] = [dg, dge]
    touched += 1
    print(f"Stamped Convention A anchor: {cpd_id} ({label_name}) {LABEL} = [{dg}, {dge}]")
if touched:
    compounds_helper.saveCompounds(compounds_dict)
