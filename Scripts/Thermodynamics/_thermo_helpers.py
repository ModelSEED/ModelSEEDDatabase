"""Reusable primitives for the thermodynamics update scripts.

Adding a new energy source typically means writing one short script that
composes these primitives, rather than copy-pasting another 60-100 line file.

Compound-energy sources need three things:
  1. a parser that returns a ``{lookup_key: [dg, dge]}`` dict
     (use ``parse_two_col_energy_table`` for ``id\\tdg\\tdge`` extracts
     or ``parse_gc_compound_table`` for MFAToolkit mol-analysis tables)
  2. a resolver ``(cpd, structure_type, structure, aliases) -> [dg, dge] | None``
     called only when a structure is present. Return ``None`` to skip
     writing for this compound (e.g. eQuilibrator's "structure unmapped"
     case). The "no structure at all" case is handled by ``run_compound_update``
     via its ``on_no_structure=`` flag — resolvers no longer branch on it.
  3. a label string used as the key under ``thermodynamics`` in the JSON
     output (e.g. ``"Group contribution"`` or ``"eQuilibrator"``).

The resolver pattern keeps each source's idiosyncratic decisions visible
in the per-source script while sharing the iteration/save scaffolding here.

Reaction-energy sources come in two flavors:
  - aggregate from compound energies for the same label
    (``iter_complete_reactions`` + ``sum_reaction_energy``)
  - direct lookup from a precomputed reaction table
    (``parse_two_col_energy_table`` + ``set_thermo``)

Output of every helper is intentionally bit-for-bit identical to the
pre-refactor scripts. See the README for the addition walkthrough.
"""
import csv
import os

# Sentinel marking "no energy available". Stored as a float so downstream
# math still works without a None check.
DEFAULT_DG = 10000000.0
DEFAULT_DG_DGE = [DEFAULT_DG, DEFAULT_DG]

# Preferred structure types in priority order. First match wins.
STRUCTURE_PREFERENCE = ('InChIKey', 'SMILE')


def _root(*parts):
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, '..', '..', *parts)


def biochem_path(*parts):
    """Path inside ``Biochemistry/`` relative to this script."""
    return _root('Biochemistry', *parts)


def thermo_path(*parts):
    """Path inside ``Biochemistry/Thermodynamics/`` relative to this script."""
    return _root('Biochemistry', 'Thermodynamics', *parts)


def structures_path(*parts):
    """Path inside ``Biochemistry/Structures/`` relative to this script."""
    return _root('Biochemistry', 'Structures', *parts)


def fmt2(value):
    """Round to 2-decimal float via the same formatter the originals used."""
    return float("{0:.2f}".format(float(value)))


def fmt_dg_dge(dg, dge):
    return [fmt2(dg), fmt2(dge)]


def _per_source_operator(rxn_entry, dg, dge, label=None):
    """Compute the per-source thermodynamic-direction operator for a single
    ``(dg, dge)`` pair via the cascade heuristic that belongs to ``label``
    (EQ heuristics for ``"eQuilibrator"``, GC heuristics otherwise).

    Lazy-imported from ``Estimate_Reaction_Reversibility`` so the helpers
    module stays free of a top-level dependency on the cascade module (and
    avoids any circular-import risk if that direction is ever reversed).
    Returns one of ``'>' | '<' | '=' | '?'``.

    This is what gets stored as the third element of each reaction-side
    ``thermodynamics[label]`` triple. Computing it at write time — rather
    than relying on a placeholder + downstream backfill — keeps each
    sublist's operator a function of THAT source's own energy, independent
    of the cascade's choice of top-level ``deltag``."""
    from Estimate_Reaction_Reversibility import reversibility_from_energy
    return reversibility_from_energy(rxn_entry, dg, dge, source=label)


def pick_structure(structures_dict, cpd, preference=STRUCTURE_PREFERENCE):
    """Return ``(structure_type, structure_string)`` honoring ``preference``;
    ``(None, None)`` if the compound has no matching structure entry."""
    if cpd not in structures_dict:
        return None, None
    for stype in preference:
        if stype in structures_dict[cpd]:
            return stype, next(iter(structures_dict[cpd][stype]))
    return None, None


def set_thermo(entry, label, value):
    """Assign ``entry['thermodynamics'][label] = value``, initializing the
    ``thermodynamics`` dict only when it is not already one. Preserves the
    pre-refactor write order so the JSON output is unchanged."""
    thermo = entry.get('thermodynamics')
    if not isinstance(thermo, dict):
        thermo = dict()
        entry['thermodynamics'] = thermo
    thermo[label] = value


def has_thermo(entry, label):
    """True iff ``entry`` carries a non-sentinel energy for ``label``."""
    if not isinstance(entry, dict):
        return False
    thermo = entry.get('thermodynamics')
    return (isinstance(thermo, dict)
            and label in thermo
            and thermo[label][0] != DEFAULT_DG)


def parse_two_col_energy_table(path):
    """Parse the standard eQuilibrator-style ``id\\tdg\\tdge[...]`` extract.

    Skips rows whose dg field is ``nan`` or contains the substring
    ``energy`` (matches the ``Unable to retrieve energy`` sentinel rows the
    Retrieve_* scripts emit). Returns ``{id: [dg, dge]}`` with 2-decimal
    formatting applied at parse time so callers can compare floats directly.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            array = line.strip().split('\t')
            if len(array) < 3:
                continue
            if 'energy' in array[1] or array[1] == 'nan':
                continue
            out[array[0]] = fmt_dg_dge(array[1], array[2])
    return out


def parse_modelseed_energy_table(path, id_col, dg_col, err_col,
                                 status_col='status', ok='ok'):
    """Parse a regenerated eQuilibrator table produced directly from ModelSEED
    structures (``ModelSEED_{Reaction,Compound}_Energies.tsv``).

    These supersede the ``MetaNetX_*_Energies.tbl`` extracts. The difference
    that matters here is the ``status`` column: the generator emits a row for
    every ModelSEED reaction/compound and says why it has no energy
    (``outside CC span``, ``unbalanced``, ``compound not in cache``,
    ``translocation only``) rather than omitting it or, as the MetaNetX-era
    retrieval did, reporting a number anyway. Only ``ok`` rows carry an energy.

    Returns ``{id: [dg, dge]}`` in kcal/mol over the ``ok`` rows only, with the
    same 2-decimal formatting :func:`parse_two_col_energy_table` applies, so
    callers can compare floats directly.

    The leading ``# cache=... params=... p_h=...`` provenance line records the
    compound cache, the component-contribution parameter set and the conditions
    the numbers were computed at. It is skipped here but must not be dropped
    from the file -- it is the only record of which cache produced these values.
    """
    out = {}
    with open(path) as fh:
        first = fh.readline()
        if not first.startswith('#'):
            fh.seek(0)
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if row.get(status_col) != ok:
                continue
            try:
                out[row[id_col]] = fmt_dg_dge(row[dg_col], row[err_col])
            except (KeyError, TypeError, ValueError):
                continue
    return out


def run_reaction_table_update(reactions_helper, label, energy_table):
    """Write ``thermodynamics[label]`` from a regenerated energy table,
    REMOVING the key where the table has no energy.

    This differs from :func:`run_reaction_lookup_update`, which leaves absent
    reactions untouched. That behaviour is right for an incremental refresh of
    the same pipeline, and wrong for a regeneration: a reaction that used to
    carry a value and no longer earns one must lose it, not keep a stale number
    from the superseded pipeline sitting under a current-looking key.

    The operator is computed here, at write time, from THIS source's own energy
    and rule set -- so the stored triple is internally consistent the moment it
    lands, and the later ``Estimate_Reaction_Reversibility.py`` pass is a
    confirmation rather than a repair."""
    reactions_dict = reactions_helper.loadReactions()
    stats = {'updated': 0, 'added': 0, 'removed': 0, 'absent': 0}
    for rxn in sorted(reactions_dict.keys()):
        entry = reactions_dict[rxn]
        fresh = energy_table.get(rxn)
        thermo = entry.get('thermodynamics')
        had = isinstance(thermo, dict) and label in thermo
        if fresh is None:
            if had:
                del thermo[label]
                stats['removed'] += 1
            else:
                stats['absent'] += 1
            continue
        dg, dge = fresh[0], fresh[1]
        op = _per_source_operator(entry, dg, dge, label)
        set_thermo(entry, label, [dg, dge, op])
        stats['updated' if had else 'added'] += 1

    for key in ('added', 'updated', 'removed', 'absent'):
        print("  %-10s %7d" % (key, stats[key]))
    print("Saving reactions")
    reactions_helper.saveReactions(reactions_dict)
    return stats


def run_compound_table_update(compounds_helper, label, energy_table):
    """Compound-side counterpart of :func:`run_reaction_table_update`.

    Compounds store ``[dg, dge]`` with no third element: a formation energy has
    no direction, so writing an operator slot would invent a field the schema
    does not have."""
    compounds_dict = compounds_helper.loadCompounds()
    stats = {'updated': 0, 'added': 0, 'removed': 0, 'absent': 0}
    for cpd in sorted(compounds_dict.keys()):
        entry = compounds_dict[cpd]
        fresh = energy_table.get(cpd)
        thermo = entry.get('thermodynamics')
        had = isinstance(thermo, dict) and label in thermo
        if fresh is None:
            if had:
                del thermo[label]
                stats['removed'] += 1
            else:
                stats['absent'] += 1
            continue
        set_thermo(entry, label, [fresh[0], fresh[1]])
        stats['updated' if had else 'added'] += 1

    for key in ('added', 'updated', 'removed', 'absent'):
        print("  %-10s %7d" % (key, stats[key]))
    print("Saving compounds")
    compounds_helper.saveCompounds(compounds_dict)
    return stats


def parse_gc_compound_table(thermo_root, sources=("KEGG", "MetaCyc"),
                            processes=("Charged", "Original")):
    """Parse Group-Contribution mol-analysis tables under ``thermo_root``.

    Source files live at
    ``{thermo_root}/ModelSEED/{source}_{process}_MolAnalysis.tbl`` and carry
    dg/dge in columns 7 and 8 (0-indexed). Entries are stored as ``str`` for
    the NoGroup-recovery comparison to work the same way as the original
    (``'10000000.00'`` vs ``'1e+07'``)."""
    table = {}
    for source in sources:
        for process in processes:
            fname = os.path.join(
                thermo_root, 'ModelSEED', f'{source}_{process}_MolAnalysis.tbl')
            with open(fname) as fh:
                for line in fh:
                    array = line.strip().split('\t')
                    if len(array) < 9:
                        continue
                    key = array[0]
                    dg_str = "{0:.2f}".format(float(array[7]))
                    dge_str = "{0:.2f}".format(float(array[8]))
                    if key not in table:
                        table[key] = {'dg': dg_str, 'dge': dge_str}
                        continue
                    # NoGroup recovery: a handful of protonated mol files get
                    # tagged 'NoGroup' by MFAToolkit and their dg sentinels
                    # out. If we previously stored the sentinel and the
                    # current row carries a real energy, prefer the real one.
                    if (table[key]['dg'] == "10000000.00"
                            and array[7] != "1e+07"):
                        table[key] = {'dg': dg_str, 'dge': dge_str}
    return table


def parse_curated_structure_aliases(all_structures_path):
    """Parse ``All_ModelSEED_Structures.txt`` into
    ``{cpd: {structure: {alias: 1}}}``.

    Used by the GC compound updater to drop alias entries that came from
    structures other than the curated one for each compound — preserving the
    provenance trim from the original script."""
    curated = {}
    with open(all_structures_path) as fh:
        for line in fh:
            array = line.strip().split('\t')
            if len(array) < 8:
                continue
            cpd, structure, alias = array[0], array[7], array[3]
            curated.setdefault(cpd, {}).setdefault(structure, {})[alias] = 1
    return curated


def parse_mnx_inchikey_index(structures_in_modelseed_and_equilibrator_path):
    """Parse the MetaNetX file mapping MNX id -> InChIKey, building a
    deprotonated (2-segment) InChIKey -> list-of-MNX-ids index. Matches
    the eQuilibrator updater's original lookup table exactly."""
    index = {}
    with open(structures_in_modelseed_and_equilibrator_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            mnx, inchikey = line.split('\t')
            inchikey = "-".join(inchikey.split('-')[0:2])
            bucket = index.setdefault(inchikey, [])
            if mnx not in bucket:
                bucket.append(mnx)
    return index


def lowest_energy_eq_style(candidate_keys, lookup):
    """List-iteration variant matching the eQuilibrator compound path.
    First strictly-lower-dg candidate wins; tie keeps the earlier one."""
    best = list(DEFAULT_DG_DGE)
    for key in candidate_keys:
        entry = lookup.get(key)
        if entry is None:
            continue
        if entry[0] < best[0]:
            best = entry
    return best


def lowest_energy_gc_style(candidate_keys, lookup):
    """Dict-bucket variant matching the GC compound path. Duplicate dg
    values collapse to the *last*-seen dge before the minimum search runs."""
    by_dg = {}
    for key in candidate_keys:
        entry = lookup.get(key)
        if entry is None:
            continue
        by_dg[float(entry['dg'])] = float(entry['dge'])
    best_dg, best_dge = DEFAULT_DG, DEFAULT_DG
    for dg, dge in by_dg.items():
        if dg < best_dg:
            best_dg, best_dge = dg, dge
    return [best_dg, best_dge]


def mean_energy_gc_style(candidate_keys, lookup):
    """Mean-averaging variant: for a set of accepted-structure aliases that
    all nominally represent the same ModelSEED compound, return the
    arithmetic mean of their ΔG values plus a variance-inflated error.

    ``variance-inflated error = sqrt(mean(σᵢ²) + var(dgᵢ))``. The first term
    is the average per-alias measurement error; the second is the spread
    across aliases (which itself signals ambiguity between the source
    mol files' raw representations of the same InChIKey). When aliases
    agree exactly, var=0 and the error reduces to the mean of the
    individual errors."""
    values = []
    errors = []
    for key in candidate_keys:
        entry = lookup.get(key)
        if entry is None:
            continue
        dg = float(entry['dg'])
        if abs(dg) >= DEFAULT_DG:
            continue
        values.append(dg)
        errors.append(float(entry['dge']))
    if not values:
        return [DEFAULT_DG, DEFAULT_DG]
    if len(values) == 1:
        return [round(values[0], 2), round(errors[0], 2)]
    mean_dg = sum(values) / len(values)
    mean_var = sum(e * e for e in errors) / len(errors)
    var_dg = sum((v - mean_dg) ** 2 for v in values) / len(values)
    combined_err = (mean_var + var_dg) ** 0.5
    return [round(mean_dg, 2), round(combined_err, 2)]


def iter_compounds(compounds_dict, structures_dict,
                   preference=STRUCTURE_PREFERENCE):
    """Yield ``(cpd, entry, structure_type, structure, aliases)`` in
    sorted-id order (the order all original update scripts use).

    ``structure_type``, ``structure`` and ``aliases`` are ``None`` when
    no structure is present. ``aliases`` is the alias list attached to
    the picked structure (a list of strings), pre-extracted so resolvers
    don't have to walk ``structures_dict`` by hand."""
    for cpd in sorted(compounds_dict.keys()):
        stype, structure = pick_structure(structures_dict, cpd, preference)
        if structure is None:
            yield cpd, compounds_dict[cpd], None, None, None
            continue
        aliases = structures_dict[cpd][stype][structure].get('alias', [])
        yield cpd, compounds_dict[cpd], stype, structure, aliases


def run_compound_update(compounds_helper, label, resolve_energy,
                        on_no_structure='default',
                        preference=STRUCTURE_PREFERENCE):
    """Shared scaffolding for compound-energy updates.

    ``resolve_energy(cpd, structure_type, structure, aliases)`` is called
    only when the compound has a structure. It returns either a
    ``[dg, dge]`` pair to write or ``None`` to skip this compound.

    ``on_no_structure`` controls what happens when the compound has no
    structure at all:
      - ``'default'`` (GC behavior): write the sentinel ``DEFAULT_DG_DGE``
      - ``'skip'``: leave the compound's thermodynamics untouched

    Sources whose structure-present resolver might also want to skip
    (eQuilibrator: structure exists but isn't in MetaNetX) just return
    ``None`` from the resolver — that path is independent of
    ``on_no_structure``."""
    compounds_dict = compounds_helper.loadCompounds()
    structures_dict = compounds_helper.loadStructures(
        ["SMILE", "InChIKey"], ["ModelSEED"])

    for cpd, entry, stype, structure, aliases in iter_compounds(
            compounds_dict, structures_dict, preference):
        if structure is None:
            if on_no_structure == 'skip':
                continue
            energy = list(DEFAULT_DG_DGE)
        else:
            energy = resolve_energy(cpd, stype, structure, aliases)
            if energy is None:
                continue
        set_thermo(entry, label, energy)

    print("Saving compounds")
    compounds_helper.saveCompounds(compounds_dict)


def eligible_compounds_for_label(compounds_dict, label):
    """Return ``{cpd: 1}`` for compounds carrying a usable energy for
    ``label``, expanding obsolete compounds via ``linked_compound``. Used
    by the per-compound-aggregating reaction updaters (GC today; any
    future source that stores per-compound energies)."""
    eligible = {}
    for cpd, entry in compounds_dict.items():
        if not has_thermo(entry, label):
            continue
        eligible[cpd] = 1
        if entry['is_obsolete']:
            for link in entry['linked_compound'].split(';'):
                if link in eligible:
                    continue
                if has_thermo(compounds_dict.get(link), label):
                    eligible[link] = 1
    return eligible


def sum_reaction_energy(stoichiometry, compounds_dict, label, rxn_id_for_warn):
    """Stoichiometric sum of per-compound energies stored under ``label``.

    Variance is summed in quadrature (matches the original). Returns a
    2-element ``[dg, dge]`` pair; callers stamp the per-source operator
    onto the result themselves (see :func:`run_reaction_aggregation_update`).
    ``rxn_id_for_warn`` is only used for the warning print when a reagent
    is unexpectedly absent."""
    dg_sum = 0.0
    dge_sq_sum = 0.0
    for rgt in stoichiometry:
        cpd = rgt['compound']
        thermo = compounds_dict.get(cpd, {}).get('thermodynamics')
        if not isinstance(thermo, dict) or label not in thermo:
            print("Warning: wrong reaction: " + rxn_id_for_warn)
            continue
        dg, dge = thermo[label][0], thermo[label][1]
        dg_sum += dg * rgt['coefficient']
        dge_sq_sum += (dge * rgt['coefficient']) ** 2
    return fmt_dg_dge(dg_sum, dge_sq_sum ** 0.5)


def run_reaction_aggregation_update(reactions_helper, compounds_helper, label):
    """Sum per-compound energies stored under ``label`` across each
    non-EMPTY reaction's stoichiometry; write a 3-element
    ``[dg, dge, operator]`` triple where ``operator`` is THIS source's own
    thermodynamic direction (computed from the same heuristic the
    canonical reversibility cascade uses, applied to this source's dG).

    Writes the default sentinel triple when any reagent lacks an energy
    under ``label`` — matching the original 'all-or-nothing' GC reaction
    semantics. ``reversibility_from_energy`` returns ``'?'`` for sentinel
    inputs, so sentinel entries naturally end up ``[SENTINEL, SENTINEL, '?']``."""
    compounds_dict = compounds_helper.loadCompounds()
    eligible = eligible_compounds_for_label(compounds_dict, label)

    reactions_dict = reactions_helper.loadReactions()
    for rxn, rxn_entry in reactions_dict.items():
        if rxn_entry['status'] == 'EMPTY':
            continue
        rgts = rxn_entry['stoichiometry']
        if all(rgt['compound'] in eligible for rgt in rgts):
            dg, dge = sum_reaction_energy(rgts, compounds_dict, label, rxn)
        else:
            dg, dge = DEFAULT_DG, DEFAULT_DG
        op = _per_source_operator(rxn_entry, dg, dge, label)
        set_thermo(rxn_entry, label, [dg, dge, op])

    print("Saving reactions")
    reactions_helper.saveReactions(reactions_dict)


def run_reaction_lookup_update(reactions_helper, label, energy_table):
    """Overwrite ``thermodynamics[label]`` from a precomputed ``{rxn: [dg, dge]}``
    table; write a 3-element ``[dg, dge, operator]`` triple where the
    operator is THIS source's own thermodynamic direction (computed via
    :func:`_per_source_operator`). Reactions absent from the table are
    left untouched — matching the original eQuilibrator reaction-updater
    behavior. The parser itself stays 2-element because it is shared with
    the compound-side resolvers, which remain 2-element."""
    reactions_dict = reactions_helper.loadReactions()
    for rxn in sorted(reactions_dict.keys()):
        if rxn not in energy_table:
            continue
        dg, dge = energy_table[rxn][0], energy_table[rxn][1]
        op = _per_source_operator(reactions_dict[rxn], dg, dge, label)
        set_thermo(reactions_dict[rxn], label, [dg, dge, op])

    print("Saving reactions")
    reactions_helper.saveReactions(reactions_dict)
