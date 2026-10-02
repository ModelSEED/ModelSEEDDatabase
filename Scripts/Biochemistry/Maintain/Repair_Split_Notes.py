#!/usr/bin/env python
"""Repair reaction notes that were split one character per element.

Adjust_Reaction_Water.py appended its marker with

    Reactions_Dict[rxn]["notes"] += "|WB"

but `notes` is a LIST, so `list += str` extended it one character at a time and
wrote ['GCP', 'EQP', '|', 'W', 'B'] where ['GCP', 'EQP', 'WB'] was meant. The
script was corrected on 2026-09-21; this repairs the records it had already
written. 15 reactions are affected.

The repair is the inverse of the bug: a '|' element followed by a run of
single-character elements is the fragmented form of one note, so the '|' is
dropped and the run is rejoined. The result is checked against the note
vocabulary rather than trusted, so a run that rejoins to something unrecognised
is reported and left alone.

Also drops empty-string notes, which carry no information and appear in the
same records.

    Repair_Split_Notes.py          dry run, report only
    Repair_Split_Notes.py save     write the repaired records
"""

if __name__ == "__main__":
    import argparse as _argparse
    _p = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter)
    _p.add_argument("modes", nargs="*", choices=("save", []),
                    help="save: write the repaired records (default is a dry run)")
    _p.parse_args()

import sys
sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions

# Every note this database uses. A rejoined run must land in here.
VALID_NOTES = {"EQP", "EQC", "EQU", "GCC", "GCP", "HB", "WB", "NB"}

dry_run = "save" not in sys.argv

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()


def repair(notes):
    """Rejoin character-split runs; drop empty strings. Returns (notes, changed)."""
    out = []
    i = 0
    while i < len(notes):
        if notes[i] == "|":
            j = i + 1
            run = []
            while j < len(notes) and len(notes[j]) == 1 and notes[j] != "|":
                run.append(notes[j])
                j += 1
            joined = "".join(run)
            if joined in VALID_NOTES:
                out.append(joined)
                i = j
                continue
            # Not a note we recognise -- leave the fragments untouched and let
            # the caller report it rather than inventing a value.
            out.append(notes[i])
            i += 1
            continue
        if notes[i] != "":
            out.append(notes[i])
        i += 1
    return out, out != list(notes)


repaired = 0
unrecognised = []
for rxn in sorted(reactions_dict.keys()):
    notes = reactions_dict[rxn].get("notes") or []
    if not isinstance(notes, list):
        continue
    new_notes, changed = repair(notes)
    if not changed:
        continue
    if "|" in new_notes:
        unrecognised.append((rxn, notes))
        continue
    print("Repairing " + rxn + ": " + str(notes) + " -> " + str(new_notes))
    reactions_dict[rxn]["notes"] = new_notes
    repaired += 1

if unrecognised:
    print("\nLeft alone, rejoined to something outside the vocabulary:")
    for rxn, notes in unrecognised:
        print("  " + rxn + ": " + str(notes))

print("\nRepaired " + str(repaired) + " reactions")
if repaired > 0:
    if dry_run:
        print("Dry run -- pass 'save' to write these changes")
    else:
        reactions_helper.saveReactions(reactions_dict)
        print("Saved")
