#!/usr/bin/env python
"""Write Group-Contribution compound energies into the BiochemPy compound JSON.

Source: ``Biochemistry/Thermodynamics/ModelSEED/ModelSEED_GroupContribution.tsv``
-- MFAToolkit's group-contribution method run directly against the structures
ModelSEED curates, keyed by ``cpd`` id.

This replaces the alias-mediated route. Previously GC decomposed KEGG and
MetaCyc mol-file corpora, keyed by *source* accession, and reached a ModelSEED
compound by matching aliases -- taking the mean across whichever source
molecules shared an alias, with the error widened by their disagreement. The
energy therefore came from whichever molecule the alias pointed at, which is
not necessarily the molecule we hold. That is the same indirection removed from
the eQuilibrator path in the 2026-08 regeneration; with it gone, all three
thermodynamic sources derive from ModelSEED structures.

Keyed by our own id there is exactly one structure per compound, so the
mean-of-aliases resolver and its variance inflation no longer have anything to
operate on and are retired. Note what that costs: for the handful of compounds
whose aliases genuinely disagreed (pteridines, flavonoids -- spreads of 1-24
kcal/mol), the old error bar advertised that disagreement. The new one does
not. The disagreement has not been resolved, it has been made invisible; what
remains is the group-contribution uncertainty for the one structure we chose.

Convention A throughout -- see the anchor table below and
Biochemistry/Thermodynamics/README.md."""

if __name__ == "__main__":
    # Validate arguments BEFORE importing anything or touching the database.
    # These scripts mutate the database, and without this an unknown flag or a
    # mistyped mode was silently ignored and the script ran with its defaults:
    # asking Estimate_Reaction_Reversibility.py for --help rewrote 122 files.
    # Placed above the imports so --help works even where a dependency is
    # missing from the path.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds
import _thermo_helpers as th

LABEL = 'Group contribution'

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
gc = th.parse_modelseed_gc_table(
    th.thermo_path('ModelSEED', 'ModelSEED_GroupContribution.tsv'))

print("%d compounds with a group-contribution energy" % len(gc))


def resolve(cpd, stype, structure, aliases):
    """Direct lookup by ModelSEED id. No alias walk, no averaging.

    A compound whose structure the group decomposer declined gets the SENTINEL,
    not a skip. Returning None would leave whatever the alias route had written
    still sitting in the record -- a value derived from some KEGG or MetaCyc
    molecule, under a label that now claims to mean "computed from the
    structure we hold". 99 compounds were in exactly that state on the first
    pass; the point of this change is that no such value survives.
    """
    return gc.get(cpd, list(th.DEFAULT_DG_DGE))


# on_no_structure='default' preserves the historical sentinel for compounds
# with no structure at all. Changing that is a separate question (the three
# sources currently disagree about how to say "no value"); this commit changes
# only where the energy comes from.
th.run_compound_update(compounds_helper, LABEL, resolve,
                       on_no_structure='default')

# Post-run: stamp all Convention A anchors. MFAToolkit cannot decompose these
# standalone -- water comes back "unlabeled atoms:3", not an energy -- so their
# ΔfG values live only in the cue database and are injected here at the
# compound level. Unchanged by the move to direct lookup: the reason these need
# injecting was never the alias route, it is that the group decomposer has
# nothing to say about them.
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
