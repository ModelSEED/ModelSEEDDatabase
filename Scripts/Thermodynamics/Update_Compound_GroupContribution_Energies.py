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
    # Argument guard -- see "The argument guard" in Scripts/README.md.
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

# Every compound MFAToolkit was RUN on, whatever the outcome. parse_modelseed_
# gc_table() returns status==ok rows only, so it cannot distinguish "declined"
# from "never submitted" -- and that distinction is the whole point here.
SUBMITTED = set()
with open(th.thermo_path('ModelSEED', 'ModelSEED_GroupContribution.tsv')) as _fh:
    for _line in _fh:
        if _line.startswith('#'):
            continue
        _id = _line.split('\t', 1)[0]
        if _id.startswith('cpd'):
            SUBMITTED.add(_id)

print("%d compounds with a group-contribution energy" % len(gc))
print("%d compounds submitted to MFAToolkit (energy or stated refusal)" % len(SUBMITTED))


def resolve(cpd, stype, structure, aliases):
    """Direct lookup by ModelSEED id. No alias walk, no averaging.

    A compound whose structure the group decomposer declined gets the SENTINEL,
    not a skip. Returning None would leave whatever the alias route had written
    still sitting in the record -- a value derived from some KEGG or MetaCyc
    molecule, under a label that now claims to mean "computed from the
    structure we hold". 99 compounds were in exactly that state on the first
    pass; the point of this change is that no such value survives.
    """
    # Three outcomes (2026-09-15):
    #   energy in the table          -> the value
    #   submitted but declined       -> SENTINEL (GC ran and could not finish)
    #   never submitted              -> None, write nothing
    if cpd in gc:
        return gc[cpd]
    if cpd in SUBMITTED:
        return list(th.DEFAULT_DG_DGE)
    return None


# on_no_structure='skip' (2026-09-15, was 'default'): a compound with no
# structure was never submitted to the group decomposer, so it gets no entry
# rather than a sentinel. This resolves the question the previous comment here
# deferred -- the three sources now agree on how to say "no value": absence.
th.run_compound_update(compounds_helper, LABEL, resolve,
                       on_no_structure='skip')

# 'skip' only declines to WRITE; it leaves any sentinel a previous run left
# behind. Purge those explicitly, sparing the Convention A anchors stamped
# below (MFAToolkit declines them, but chemistry supplies their values).
_anchors = {a[0] for a in CONVENTION_A_ANCHORS}
_cd = compounds_helper.loadCompounds()
_purged = 0
for _cid, _entry in _cd.items():
    if _cid in SUBMITTED or _cid in _anchors:
        continue
    if th.drop_thermo(_entry, LABEL):
        _purged += 1
if _purged:
    print(f"Purged {_purged:,} GC entries for compounds never submitted to MFAToolkit")
    compounds_helper.saveCompounds(_cd)

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
