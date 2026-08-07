#!/usr/bin/env python
"""Estimate reaction reversibility (``>``, ``<``, ``=``, or ``?``) from the
stored thermodynamic energies and write it back into the reactions JSON.

The cascade is now **composable**: heuristics and energy sources live in
``reversibility_heuristics`` as plug-in pieces, and this module wires the
historical default rule set (``DEFAULT_HEURISTICS``) to the top-level energy
source. To use a different rule set, build your own heuristic list and/or
energy source and call ``run_reversibility`` directly, e.g.::

    from reversibility_heuristics import (
        run_reversibility, DEFAULT_HEURISTICS, per_source_energy,
        make_ln_reversibility_index_heuristic)
    rules = DEFAULT_HEURISTICS[:-1] + [make_ln_reversibility_index_heuristic(ln_ri)]
    status, op, label = run_reversibility(rxn_entry, per_source_energy("eQuilibrator"), rules)

The default ``estimate_one`` / ``reversibility_from_energy`` behaviour (and the
generated reports) are byte-for-byte unchanged.

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
)


# ---------------------------------------------------------------------------
# Cascade entry points (thin wrappers over the composable core)
# ---------------------------------------------------------------------------
def _cascade(rxn_entry, rxn_dg, rxn_dge):
    """Run the default heuristic cascade against an explicit ``(rxn_dg, rxn_dge)``
    pair and return ``(status_label, operator)``. Kept for callers that import it
    directly; equivalent to ``run_reversibility`` with ``explicit_energy`` and
    ``DEFAULT_HEURISTICS``."""
    status, operator, _ = run_reversibility(
        rxn_entry, explicit_energy(rxn_dg, rxn_dge), DEFAULT_HEURISTICS)
    return status, operator


def estimate_one(rxn_entry, db_level):
    """Returns ``(status_label, thermoreversibility, source_label)`` for one
    reaction, using the top-level energy source + the default cascade.

    ``source_label`` is the Thermodynamics subkey whose energy fed the estimate
    (or ``None`` for empty/incomplete, or when the unfiltered run's top-level
    energy did not match a sublist exactly)."""
    if rxn_entry['status'] == "EMPTY":
        return "Empty", "?", None

    status, thermoreversibility, source_label = run_reversibility(
        rxn_entry, top_level_energy(db_level), DEFAULT_HEURISTICS)
    if status is None:  # no usable energy -> incomplete fallback
        status, thermoreversibility = _incomplete_decision(rxn_entry, db_level)
        return status, thermoreversibility, None
    return status, thermoreversibility, source_label


def reversibility_from_energy(rxn_entry, rxn_dg, rxn_dge):
    """Compute the thermodynamic direction operator for a single per-source
    ``(dg, dge)`` pair without the source-eligibility filter or the top-level
    deltag pick. Returns one of ``'>'`` / ``'<'`` / ``'='`` / ``'?'``.

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

    _status, operator = _cascade(rxn_entry, dg, dge)
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
    if len(argv) > 1 and argv[1] in ('EQ', 'GC', 'DGP'):
        return argv[1]
    return ''


def main():
    db_level = _parse_db_level(sys.argv)
    helper = Reactions()
    reactions_dict = helper.loadReactions()

    report = {}
    for rxn in sorted(reactions_dict.keys()):
        rxn_entry = reactions_dict[rxn]
        # The cascade-winner source label is intentionally ignored here:
        # per-source operators are written at energy-table time by
        # ``_thermo_helpers`` (each using THAT source's own dG). This step
        # only updates the canonical top-level reversibility.
        status, thermoreversibility, _ = estimate_one(rxn_entry, db_level)
        report[rxn] = [status, rxn_entry["reversibility"], thermoreversibility]
        rxn_entry['reversibility'] = thermoreversibility

    _write_report(db_level, report)
    print("Saving reactions")
    helper.saveReactions(reactions_dict)


if __name__ == "__main__":
    main()
