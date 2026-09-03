#!/usr/bin/env python
"""Encoding for the per-source ``pkas`` field on a compound record.

Atom indices are deliberately absent here.

Nothing in the energy path consumes one. eQuilibrator stores
``dissociation_constants`` as a positional list of floats and only ever takes
its length, compares entries to the pH and sums a slice; equilibrator-assets
takes ``specified_pkas`` as SMILES -> list of floats; our cache builder emits
``";".join(f"{v:.4f}")``; ``resolve_pkas`` returns ``List[float]``; dGPredictor
does not use pKas at all; Group Contribution does not reference them. The only
reader of the atom slot was ``Build_Biochemistry_Object.py``, and 64.2% of the
keys it published pointed at carbon, phosphorus, or an atom index beyond the
molecule -- because Marvin reorders atoms on import, so its indices never
described our structures.

The indices are not discarded, only relocated. They stay in the per-source
files under ``Biochemistry/Structures/<source>/pkas/``, each in the atom space
it was actually computed in: Marvin's under KEGG / MetaCyc / ChEBI / Rhea keyed
by external id, MolGpKa's under ModelSEED keyed by cpd id. That is the
provenance record, and it is the only place an index means anything.

Encoding is ``<fragment>:<value>``. The fragment index is kept because it is
real -- a salt contributes pKas from its counter-ion, and a consumer wanting
the parent molecule needs to filter them out. Ordering of macroscopic values is
the order of the tokens.

Each source entry carries ``kind``:

  microscopic   per-site values predicted on one protonation state, as Marvin
                and MolGpKa report. Not a ladder; slicing or summing them is
                meaningless.
  macroscopic   sequential dissociation constants of the whole molecule, as
                measured values are. This is what eQuilibrator's transform
                actually requires.

Keeping that distinction explicit is the point: eQuilibrator reads any pKa list
as a macroscopic ladder, so handing it microscopic values is a category error,
not merely an inaccuracy. It moved one LPS reaction by 355 kcal/mol.
"""

MICROSCOPIC = "microscopic"
MACROSCOPIC = "macroscopic"


def encode(pkas, fragment=1):
    """pkas: iterable of values. For a macroscopic set, in dissociation order."""
    return ";".join(f"{fragment}:{float(p):.2f}" for p in pkas)


def decode(text, fragment=1):
    """Return the list of pKa values for one fragment, in stored order."""
    out = []
    for tok in (text or "").split(";"):
        if not tok:
            continue
        parts = tok.split(":")
        if len(parts) != 2:
            raise ValueError(
                f"malformed pKa token {tok!r}; expected <fragment>:<value>. "
                f"Three-field tokens are the retired atom-indexed form and "
                f"live in Biochemistry/Structures/<source>/pkas/ instead."
            )
        if int(parts[0]) == fragment:
            out.append(float(parts[1]))
    return out


def ladder(entry, fragment=1):
    """Macroscopic pKas as an ordered list -- the form eQuilibrator wants.

    ``entry`` is one source dict, e.g. ``{"kind": ..., "pKa": ..., "pKb": ...}``.
    Raises on a microscopic entry: silently treating per-site values as
    sequential dissociations is the bug this module exists to prevent.
    """
    if entry.get("kind") != MACROSCOPIC:
        raise ValueError(
            f"refusing to build a ladder from a {entry.get('kind')!r} pKa set: "
            f"these are per-site values on one protonation state, not "
            f"sequential dissociations of the molecule"
        )
    return decode(entry.get("pKa"), fragment) + decode(entry.get("pKb"), fragment)


def migrate(text, fragment=1):
    """Convert a stored three-field ``<fragment>:<atom>:<value>`` string to the
    two-field form, dropping the atom index."""
    out = []
    for tok in (text or "").split(";"):
        if not tok:
            continue
        parts = tok.split(":")
        if len(parts) == 2:
            out.append(tok)
        elif len(parts) == 3:
            if int(parts[0]) == fragment:
                out.append(f"{parts[0]}:{float(parts[2]):.2f}")
        else:
            raise ValueError(f"malformed pKa token {tok!r}")
    return ";".join(out)


if __name__ == "__main__":
    lit = {"kind": MACROSCOPIC, "pKa": encode([2.15, 7.20, 12.35]), "pKb": ""}
    print("literature  :", lit["pKa"])
    print("  ladder    :", ladder(lit))

    mol = {"kind": MICROSCOPIC, "pKa": encode([2.11, 2.11, 2.11]), "pKb": ""}
    print("MolGpKa     :", mol["pKa"])
    try:
        ladder(mol)
    except ValueError as e:
        print("  ladder    : refused ->", str(e)[:58] + "...")

    old = "1:2:12.90;1:3:1.80;1:4:6.95"
    print("migrate     :", old, "->", migrate(old))
    assert decode(encode([2.15, 7.20])) == [2.15, 7.20]
    assert migrate("1:2:12.90;2:9:3.00") == "1:12.90"
    print("\nassertions passed")
