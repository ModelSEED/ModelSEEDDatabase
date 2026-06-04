#!/usr/bin/env python
"""
Populate compound['properties'] from per-source Marvin TSVs.

Each qualifying compound gains:

    "properties": {
        "Marvin": [pkas_string, pkbs_string, charge_at_ph7_int]
    }

Source cascade: KEGG > MetaCyc > ChEBI > Rhea. First source with any
Marvin data for the compound (via its aliases) wins the whole tuple,
matching how Scripts/Structures/Update_Compound_pKas.py picks the
top-level pka/pkb fields.

When two sources both have data for the same MS compound and they
disagree (semantically — see compare_field below), the disagreement
is logged to Biochemistry/Structures/_reports/Properties_Conflicts.md.
"""

import os, sys
from datetime import date
sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

PKA_DBS    = ["KEGG", "MetaCyc", "ChEBI", "Rhea"]
TOOL_LABEL = "Marvin"
REPORT_PATH = "../../Biochemistry/Structures/_reports/Properties_Conflicts.md"


def parse_logks(value_str):
    """coeff:atom:logK[;coeff:atom:logK;...] -> sorted tuple of rounded logKs.

    Atom indices are dropped because different source libraries number the
    same molecule's atoms differently; the chemistry (logK values) is what
    matters. Two-decimal rounding swallows trailing-precision noise.
    """
    if not value_str:
        return ()
    out = []
    for entry in value_str.split(';'):
        parts = entry.split(':')
        if len(parts) < 3:
            continue
        try:
            out.append(round(float(parts[2]), 2))
        except ValueError:
            continue
    return tuple(sorted(out))


def compare_field(name, a, b):
    """Return True iff two values for the same field disagree semantically."""
    if name in ('pka', 'pkb'):
        return parse_logks(a) != parse_logks(b)
    # charge: direct equality (None == None is fine; only flag a real diff)
    return a != b


def write_report(path, conflicts, n_total, n_with_props):
    by_field = {'pka': 0, 'pkb': 0, 'charge': 0}
    conflicting_cpds = set()
    for c in conflicts:
        conflicting_cpds.add(c['cpd'])
        for f in c['fields']:
            by_field[f] = by_field.get(f, 0) + 1

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write("# Compound Properties — per-source conflicts (Marvin)\n\n")
        fh.write("Cascade: KEGG > MetaCyc > ChEBI > Rhea. "
                 "The \"Picked\" column reports the value written to "
                 "`properties[\"" + TOOL_LABEL + "\"]`.\n\n")
        fh.write("Generated: " + date.today().isoformat() +
                 " by Scripts/Structures/Update_Compound_Properties.py\n\n")
        fh.write("## Summary\n\n")
        fh.write("- Total compounds scanned: " + str(n_total) + "\n")
        fh.write("- Compounds with properties populated: " + str(n_with_props) + "\n")
        fh.write("- Compounds with at least one cross-source conflict: "
                 + str(len(conflicting_cpds)) + "\n")
        fh.write("- Conflict rows by field: pKa: " + str(by_field['pka'])
                 + ", pKb: " + str(by_field['pkb'])
                 + ", charge: " + str(by_field['charge']) + "\n\n")

        if not conflicts:
            fh.write("## Conflicts\n\n_No cross-source conflicts detected._\n")
            return

        fh.write("## Conflicts\n\n")
        fh.write("One row per (compound, conflicting-source, field) tuple.\n\n")
        fh.write("| MS compound | Field | Picked value | Picked source | Conflicting source | Conflicting value |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for c in conflicts:
            for f in c['fields']:
                pv = c['picked'].get(f)
                ov = c['other'].get(f)
                fh.write("| " + c['cpd']
                         + " | " + f
                         + " | " + ('' if pv is None else str(pv))
                         + " | " + c['picked_db']
                         + " | " + c['other_db']
                         + " | " + ('' if ov is None else str(ov))
                         + " |\n")


def main():
    compounds_helper = Compounds()
    compounds_dict   = compounds_helper.loadCompounds()
    structures_dict  = compounds_helper.loadStructures(
        ["SMILE", "InChI", "InChIKey"], ["ModelSEED"])
    aliases_dict     = compounds_helper.loadMSAliases()

    per_source_pkas    = compounds_helper.loadPerSourcePkas(PKA_DBS)
    per_source_charges = compounds_helper.loadPerSourceCharges(PKA_DBS, 7)

    # Strip existing properties field so re-runs are idempotent
    for cpd in compounds_dict:
        compounds_dict[cpd].pop('properties', None)

    conflicts    = []
    n_with_props = 0

    for cpd in structures_dict:
        per_source = {}    # db -> {'pka':str, 'pkb':str, 'charge':int|None}
        for db in PKA_DBS:
            if db not in aliases_dict.get(cpd, {}):
                continue
            for alias in aliases_dict[cpd][db]:
                pka_entry = per_source_pkas.get((db, alias), {})
                charge    = per_source_charges.get((db, alias))
                if pka_entry or charge is not None:
                    per_source[db] = {
                        'pka':    pka_entry.get('pKa', ''),
                        'pkb':    pka_entry.get('pKb', ''),
                        'charge': charge,
                    }
                    break   # one alias per db is enough — same chemistry

        if not per_source:
            continue

        picked_db = next(db for db in PKA_DBS if db in per_source)
        picked    = per_source[picked_db]

        for db in PKA_DBS:
            if db == picked_db or db not in per_source:
                continue
            other = per_source[db]
            diffs = [f for f in ('pka', 'pkb', 'charge')
                     if compare_field(f, picked[f], other[f])]
            if diffs:
                conflicts.append({
                    'cpd':       cpd,
                    'picked_db': picked_db,
                    'other_db':  db,
                    'fields':    diffs,
                    'picked':    picked,
                    'other':     other,
                })

        compounds_dict[cpd]['properties'] = {
            TOOL_LABEL: [
                picked['pka'],
                picked['pkb'],
                picked['charge'] if picked['charge'] is not None else '',
            ],
        }
        n_with_props += 1

    print("Saving compounds")
    compounds_helper.saveCompounds(compounds_dict)

    print("Writing conflict report to " + REPORT_PATH)
    write_report(REPORT_PATH, conflicts, len(compounds_dict), n_with_props)

    print("Done. " + str(n_with_props) + " compounds populated; "
          + str(len({c['cpd'] for c in conflicts})) + " with conflicts.")


if __name__ == '__main__':
    main()
