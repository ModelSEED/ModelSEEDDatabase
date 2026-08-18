#!/usr/bin/env python
"""Estimate reaction reversibility (``>``, ``<``, ``=``, or ``?``) from the
stored thermodynamic energies and write it back into the reactions JSON.

The cascade is **composable and source-specific**: heuristics and energy sources
live in ``reversibility_heuristics`` as plug-in pieces, and the rule set is
chosen per thermodynamic data source, because the sources do not fail the same
way.

    ./Estimate_Reaction_Reversibility.py            # top-level deltag, GC rules
    ./Estimate_Reaction_Reversibility.py GC         # Group contribution, GC rules
    ./Estimate_Reaction_Reversibility.py EQ         # eQuilibrator, EQ rules
    ./Estimate_Reaction_Reversibility.py EQ --heuristics EQ2   # eQuilibrator 2.0 rules
    ./Estimate_Reaction_Reversibility.py EQ --heuristics GC    # old behaviour

``GC`` is the default rule set: it is what every level other than ``EQ``
selects, and what any unrecognised source falls back to. The GC cascade itself
is unchanged, so ``GC`` and unfiltered runs still reproduce the historical
report byte-for-byte.

``EQ`` selects the eQuilibrator rule set (Beber 2022 uncertainty handling over
the Noor 2012 / Flamholz 2012 reversibility index) *and* switches the energy
source to ``thermodynamics['eQuilibrator']`` — see
``reversibility_heuristics.energy_source_for_level`` for why the top-level
``deltag`` was the wrong input here.

To assemble something else, call ``run_reversibility`` directly::

    from reversibility_heuristics import (
        run_reversibility, get_heuristics, per_source_energy)
    status, op, label = run_reversibility(
        rxn_entry, per_source_energy("eQuilibrator"), get_heuristics("EQ"))

The per-source ``GCC``/``EQU`` notes are no longer consulted; ``GC`` and ``EQ``
runs read directly from ``thermodynamics['Group contribution']`` and
``thermodynamics['eQuilibrator']``. After estimation, the computed direction is
appended to whichever Thermodynamics sublist supplied the energy."""
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions

# Composable cascade core. Re-imported here so the historical public surface of
# this module (constants + ``_*`` building blocks + ``estimate_one`` /
# ``reversibility_from_energy``) keeps working for existing importers
# (Update_Reaction_dGPredictor_Energies.py, _thermo_helpers, the tests).
from reversibility_heuristics import (
    # constants
    TEMPERATURE, GAS_CONSTANT, RT_CONST, FARADAY, SENTINEL_DG,
    CELL_MAX, CELL_MIN, CELL_CONC, PROTON, WATER, CO2, PROTON_WATER,
    LOW_LOCAL_CONC, ATPS_REAGENTS, ATP, PHOSPHATE_IDS, LOW_ENERGY_CPDS,
    DB_LEVEL_LABEL, DB_LEVEL_NOTE, DB_LEVEL_PRIORITY,
    # building-block helpers (re-exported for back-compat importers)
    _thermo_pair, _is_source_eligible, _energy_for, _has_gc_data,
    _incomplete_decision, _walk_stoichiometry, _stored_bounds,
    _is_atp_synthase, _abc_transporter_decision, _low_energy_points,
    # composable framework
    Context, run_reversibility, DEFAULT_HEURISTICS,
    top_level_energy, per_source_energy, explicit_energy,
    stored_bounds_heuristic, atp_synthase_heuristic, abc_transporter_heuristic,
    mmdeltag_band_heuristic, low_energy_heuristic, default_heuristic,
    make_ln_reversibility_index_heuristic,
    # source-specific rule sets
    GC_HEURISTICS, EQ_HEURISTICS, EQ2_HEURISTICS, HEURISTIC_SETS,
    DEFAULT_HEURISTIC_SET, get_heuristics, heuristics_for_source,
    energy_source_for_level,
)


# ---------------------------------------------------------------------------
# Cascade entry points (thin wrappers over the composable core)
# ---------------------------------------------------------------------------
def _cascade(rxn_entry, rxn_dg, rxn_dge, heuristics=None):
    """Run a heuristic cascade against an explicit ``(rxn_dg, rxn_dge)`` pair and
    return ``(status_label, operator)``. Kept for callers that import it
    directly; equivalent to ``run_reversibility`` with ``explicit_energy``.
    ``heuristics`` defaults to the GC rule set."""
    status, operator, _ = run_reversibility(
        rxn_entry, explicit_energy(rxn_dg, rxn_dge), heuristics or GC_HEURISTICS)
    return status, operator


def estimate_one(rxn_entry, db_level, heuristics=None, energy_source=None):
    """Returns ``(status_label, thermoreversibility, source_label)`` for one
    reaction.

    ``heuristics`` defaults to the rule set that matches ``db_level`` (GC for
    everything except ``EQ``), and ``energy_source`` to the energy that rule set
    expects. Pass either explicitly to override.

    ``source_label`` is the Thermodynamics subkey whose energy fed the estimate
    (or ``None`` for empty/incomplete, or when the unfiltered run's top-level
    energy did not match a sublist exactly)."""
    if rxn_entry['status'] == "EMPTY":
        return "Empty", "?", None

    if heuristics is None:
        heuristics = get_heuristics(db_level)   # '' / 'DGP' / unknown -> GC
    if energy_source is None:
        energy_source = energy_source_for_level(db_level)

    status, thermoreversibility, source_label = run_reversibility(
        rxn_entry, energy_source, heuristics)
    if status is None:  # no usable energy -> incomplete fallback
        status, thermoreversibility = _incomplete_decision(rxn_entry, db_level)
        return status, thermoreversibility, None
    return status, thermoreversibility, source_label


def reversibility_from_energy(rxn_entry, rxn_dg, rxn_dge, source=None):
    """Compute the thermodynamic direction operator for a single per-source
    ``(dg, dge)`` pair without the source-eligibility filter or the top-level
    deltag pick. Returns one of ``'>'`` / ``'<'`` / ``'='`` / ``'?'``.

    ``source`` is the ``thermodynamics`` subkey the pair came from (e.g.
    ``"eQuilibrator"``); it selects the rule set, defaulting to GC for every
    source without one of its own. Callers that iterate a reaction's
    ``thermodynamics`` dict should pass it — otherwise an eQuilibrator energy
    gets scored with Group-Contribution rules.

    Used by the per-source updaters (``Update_Reaction_dGPredictor_Energies.py``)
    and the operator backfill (``Add_Reaction_Thermodynamics_Operators.py``).
    Input coercion mirrors the upstream per-source updater:
      * ``rxn_entry['status'] == 'EMPTY'`` -> ``'?'``
      * ``rxn_dg`` that cannot be ``float()``-coerced (``None``, bools, NaN) -> ``'?'``
      * ``rxn_dg == SENTINEL_DG`` -> ``'?'``
      * ``rxn_dge`` that cannot be coerced -> treated as ``0.0``"""
    if isinstance(rxn_entry, dict) and rxn_entry.get('status') == 'EMPTY':
        return '?'

    if isinstance(rxn_dg, bool) or rxn_dg is None:
        return '?'
    try:
        dg = float(rxn_dg)
    except (TypeError, ValueError):
        return '?'
    if dg != dg:  # NaN
        return '?'
    if dg == SENTINEL_DG:
        return '?'

    if isinstance(rxn_dge, bool) or rxn_dge is None:
        dge = 0.0
    else:
        try:
            dge = float(rxn_dge)
        except (TypeError, ValueError):
            dge = 0.0
        if dge != dge:  # NaN
            dge = 0.0

    _status, operator = _cascade(rxn_entry, dg, dge, heuristics_for_source(source))
    return operator


# ---------------------------------------------------------------------------
# Report writer
# ---------------------------------------------------------------------------
def _write_report(db_level, report):
    """Format matches the original: GC runs drop the original-reversibility
    column from the report, EQ and unfiltered runs keep it."""
    name = "Estimated_Reaction_Reversibility_Report"
    if db_level:
        name += "_" + db_level
    name += ".txt"
    with open(name, "w") as fh:
        for rxn in sorted(report):
            row = list(report[rxn])
            if db_level == "GC":
                del row[1]
            fh.write(rxn + "\t" + "\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def _parse_db_level(argv):
    for arg in argv[1:]:
        if arg in ('EQ', 'GC', 'DGP'):
            return arg
    return ''


def _parse_heuristics(argv):
    """``--heuristics NAME`` / ``--heuristics=NAME`` override, or ``None`` to
    let the db_level pick. Rejects unknown names rather than silently
    falling back to GC, which would be easy to miss in a pipeline log."""
    for index, arg in enumerate(argv[1:], start=1):
        name = None
        if arg == '--heuristics' and index + 1 < len(argv):
            name = argv[index + 1]
        elif arg.startswith('--heuristics='):
            name = arg.split('=', 1)[1]
        if name is None:
            continue
        if name not in HEURISTIC_SETS:
            sys.exit("ERROR: unknown heuristic set %r; choose from %s"
                     % (name, ', '.join(sorted(HEURISTIC_SETS))))
        return name
    return None


def main():
    db_level = _parse_db_level(sys.argv)
    heuristics_name = _parse_heuristics(sys.argv)
    heuristics = get_heuristics(heuristics_name) if heuristics_name else None

    effective = heuristics_name or (db_level if db_level in HEURISTIC_SETS
                                    else DEFAULT_HEURISTIC_SET)
    print("Energy source: %s | heuristics: %s"
          % (db_level or 'top-level deltag', effective))

    helper = Reactions()
    reactions_dict = helper.loadReactions()

    report = {}
    for rxn in sorted(reactions_dict.keys()):
        rxn_entry = reactions_dict[rxn]
        # The cascade-winner source label is intentionally ignored here:
        # per-source operators are written at energy-table time by
        # ``_thermo_helpers`` (each using THAT source's own dG). This step
        # only updates the canonical top-level reversibility.
        status, thermoreversibility, _ = estimate_one(
            rxn_entry, db_level, heuristics=heuristics)
        report[rxn] = [status, rxn_entry["reversibility"], thermoreversibility]
        rxn_entry['reversibility'] = thermoreversibility

    _write_report(db_level, report)
    print("Saving reactions")
    helper.saveReactions(reactions_dict)


if __name__ == "__main__":
    main()
