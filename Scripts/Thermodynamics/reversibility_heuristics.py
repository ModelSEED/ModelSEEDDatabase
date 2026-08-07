"""Composable reaction-reversibility cascade: pluggable heuristics + energy sources.

The cascade is (energy_source, [heuristics]). The energy source produces
``(dg, dge, source_label)``; each heuristic is a ``(ctx: Context) -> (status, op)
| None`` callable, and the first non-None result wins.

Extend by:
  - **New energy source**: write ``(rxn_entry) -> (dg, dge, label)``.
    See :func:`top_level_energy` / :func:`per_source_energy` / :func:`explicit_energy`.
  - **New heuristic**: write ``(ctx) -> (status, op) | None`` and insert into
    the list passed to :func:`run_reversibility`, e.g.
    ``DEFAULT_HEURISTICS + [make_ln_reversibility_index_heuristic(ln_ri)]``.

``DEFAULT_HEURISTICS`` reproduces the historical fixed cascade byte-for-byte
(same order, same status strings). The regression test byte-compares the
generated reversibility report against upstream/dev, so any change here needs
to preserve those exact status strings and operator returns.
"""
from math import log

# --- Constants ------------------------------------------------------------
TEMPERATURE = 298.15
GAS_CONSTANT = 0.0019858775
RT_CONST = TEMPERATURE * GAS_CONSTANT
FARADAY = 0.023061
SENTINEL_DG = 10000000

CELL_MAX, CELL_MIN, CELL_CONC = 0.02, 0.00001, 0.001

PROTON, WATER, CO2, ATP = "cpd00067", "cpd00001", "cpd00011", "cpd00002"
PROTON_WATER = frozenset((PROTON, WATER))
LOW_LOCAL_CONC = frozenset(("cpd00007", "cpd11640"))                              # O2, H2
ATPS_REAGENTS = frozenset(("cpd00002", "cpd00008", "cpd00009",
                           "cpd00001", "cpd00067"))
PHOSPHATE_IDS = ("cpd00002", "cpd00008", "cpd00018", "cpd00009", "cpd00012")     # ATP, ADP, AMP, Pi, PPi
LOW_ENERGY_CPDS = ("cpd00011", "cpd00013", "cpd11493", "cpd00009", "cpd00012",   # CO2, NH3, ACP, Pi, PPi
                   "cpd00010", "cpd00449", "cpd00242")                            # CoA, Dihydrolipoamide, HCO3

DB_LEVEL_LABEL = {"GC": "Group contribution", "EQ": "eQuilibrator", "DGP": "dGPredictor"}
DB_LEVEL_NOTE = {"GC": "GCC", "EQ": "EQU"}
DB_LEVEL_PRIORITY = ("EQ", "GC", "DGP")

LN_RI_THRESHOLD = 6.9077552789821    # ln(1000); Noor 2012 default


# --- Energy / eligibility -------------------------------------------------
def _thermo_pair(rxn_entry, label):
    """Return ``[dg, dge]`` from ``thermodynamics[label]`` when present and
    non-sentinel, else ``None``. Length-3 entries (operator appended) tolerated."""
    thermo = rxn_entry.get('thermodynamics')
    if not isinstance(thermo, dict):
        return None
    pair = thermo.get(label)
    if not pair or pair[0] is None:
        return None
    dg = float(pair[0])
    if dg == SENTINEL_DG:
        return None
    return [dg, float(pair[1])]


def _is_source_eligible(rxn_entry, level):
    """Eligible if the structured sublist has a non-sentinel pair OR the
    legacy ``DB_LEVEL_NOTE`` flag is present in ``rxn_entry['notes']``."""
    if _thermo_pair(rxn_entry, DB_LEVEL_LABEL[level]) is not None:
        return True
    note = DB_LEVEL_NOTE.get(level)
    return note is not None and note in rxn_entry["notes"]


def _has_gc_data(rxn_entry):
    return _is_source_eligible(rxn_entry, "GC")


def _energy_for(rxn_entry, db_level):
    """Historical top-level source. Values from ``deltag``/``deltagerr``;
    ``db_level`` gates eligibility. Returns ``(dg, dge, append_label)`` or
    ``(None, None, None)``. ``append_label`` is the thermodynamics subkey
    that supplied the matching per-source pair (or None)."""
    rxn_dg = rxn_entry['deltag']
    if rxn_dg is None:
        return None, None, None
    rxn_dg = float(rxn_dg)
    if rxn_dg == SENTINEL_DG:
        return None, None, None
    rxn_dge = rxn_entry['deltagerr']
    if rxn_dge is not None:
        rxn_dge = float(rxn_dge)

    if db_level:
        if not _is_source_eligible(rxn_entry, db_level):
            return None, None, None
        label = DB_LEVEL_LABEL[db_level]
        append = label if _thermo_pair(rxn_entry, label) is not None else None
        return rxn_dg, rxn_dge, append

    for level in DB_LEVEL_PRIORITY:
        label = DB_LEVEL_LABEL[level]
        pair = _thermo_pair(rxn_entry, label)
        if pair is not None and abs(pair[0] - rxn_dg) < 1e-9:
            return rxn_dg, rxn_dge, label
    return rxn_dg, rxn_dge, None


def _incomplete_decision(rxn_entry, db_level):
    """No-usable-energy fallback. EQ inherits any prior GC reversibility."""
    if db_level == "EQ" and _has_gc_data(rxn_entry):
        return "Incomplete (GCC)", rxn_entry["reversibility"]
    return "Incomplete", "?"


# --- Stoichiometry walk (single pass, H2+H3 fixed) ------------------------
def _walk_stoichiometry(stoichiometry):
    """One pass producing every per-reaction accumulator the heuristics need."""
    rct_min = rct_max = pdt_min = pdt_max = rgt_sum = 0.0
    proton_cpts, phosphates = {}, {}
    for rgt in stoichiometry:
        cpd = rgt['compound']
        coeff = float(rgt['coefficient'])
        if cpd == PROTON:
            proton_cpts[rgt['compartment']] = 1
        if cpd in PHOSPHATE_IDS:
            phosphates[cpd] = phosphates.get(cpd, 0.0) + coeff
        if cpd in PROTON_WATER:
            continue
        if coeff < 0:
            rct_min += coeff * log(CELL_MIN)
            rct_max += coeff * log(CELL_MAX)
        else:
            pdt_min += coeff * log(CELL_MIN)
            pdt_max += coeff * log(CELL_MAX)
        local = CELL_CONC
        if cpd == CO2:
            local = 0.0001
        elif cpd in LOW_LOCAL_CONC:
            local = 0.000001
        rgt_sum += coeff * log(local)
    return {'rct_min': rct_min, 'rct_max': rct_max,
            'pdt_min': pdt_min, 'pdt_max': pdt_max,
            'rgt_sum': rgt_sum,
            'proton_cpts': proton_cpts, 'phosphates': phosphates}


def _stored_bounds(dg, dge, t):
    """(max, min) stored deltaG including concentration-range terms."""
    hi = dg + dge + RT_CONST * (t['pdt_max'] + t['rct_min'])
    lo = dg - dge + RT_CONST * (t['pdt_min'] + t['rct_max'])
    return hi, lo


def _is_atp_synthase(rxn_entry, proton_cpts):
    """Transport, multiple proton compartments, exactly the five ATPS reagents,
    only protons crossing the membrane."""
    if rxn_entry['is_transport'] != 1 or len(proton_cpts) <= 1:
        return False
    cpds_cpts = {}
    for rgt in rxn_entry['stoichiometry']:
        cpds_cpts.setdefault(rgt['compound'], []).append(rgt['compartment'])
    if len(cpds_cpts) != 5:
        return False
    for cpd, cpts in cpds_cpts.items():
        if cpd not in ATPS_REAGENTS or (len(cpts) == 2 and cpd != PROTON):
            return False
    return True


def _abc_transporter_decision(rxn_entry, phosphates):
    """Transport with ATP: direction follows the sign of the ATP coefficient."""
    if rxn_entry['is_transport'] != 1 or ATP not in phosphates:
        return None
    coeff = phosphates[ATP]
    rev = ">" if coeff < 0 else ("<" if coeff > 0 else "=")
    return f"ABCT: {coeff}", rev


def _low_energy_points(stoichiometry, phosphates):
    """Phosphate spread + low-energy-compound coefficients."""
    points = 0.0
    if ATP in phosphates and len(phosphates) > 2:
        points -= abs(min(phosphates.values()))
    for rgt in stoichiometry:
        if rgt['compound'] in LOW_ENERGY_CPDS:
            points -= float(rgt['coefficient'])
    return points


# --- Context passed to every heuristic ------------------------------------
class Context:
    """Bundled per-reaction state. ``terms`` and ``mMdeltaG`` are cached lazily
    so a chain of N heuristics costs only one stoichiometry walk."""
    __slots__ = ('rxn_entry', 'dg', 'dge', '_terms', '_mMdeltaG')

    def __init__(self, rxn_entry, dg, dge):
        self.rxn_entry, self.dg, self.dge = rxn_entry, dg, dge
        self._terms = None
        self._mMdeltaG = None

    @property
    def terms(self):
        if self._terms is None:
            self._terms = _walk_stoichiometry(self.rxn_entry['stoichiometry'])
        return self._terms

    @property
    def mMdeltaG(self):
        if self._mMdeltaG is None:
            self._mMdeltaG = self.dg + RT_CONST * self.terms['rgt_sum']
        return self._mMdeltaG


# --- Default heuristics ---------------------------------------------------
# Signature: (ctx: Context) -> (status_label, operator) | None
# First non-None wins. Extend by defining another and appending to a rules list.

def stored_bounds_heuristic(ctx):
    """MdeltaG bounds over the concentration range."""
    hi, lo = _stored_bounds(ctx.dg, ctx.dge, ctx.terms)
    if hi < 0:
        return f"MdeltaG(Max): {hi:.2f}", ">"
    if lo > 0:
        return f"MdeltaG(Min): {lo:.2f}", "<"
    return None


def atp_synthase_heuristic(ctx):
    if _is_atp_synthase(ctx.rxn_entry, ctx.terms['proton_cpts']):
        return "ATPS", "="
    return None


def abc_transporter_heuristic(ctx):
    return _abc_transporter_decision(ctx.rxn_entry, ctx.terms['phosphates'])


def mmdeltag_band_heuristic(ctx):
    if -2.0 <= ctx.mMdeltaG <= 2.0:
        return f"mMdeltaG: {ctx.mMdeltaG:.2f}", "="
    return None


def low_energy_heuristic(ctx):
    m = ctx.mMdeltaG
    pts = _low_energy_points(ctx.rxn_entry['stoichiometry'], ctx.terms['phosphates'])
    if pts * m > 2:
        return f"lowE: {m:.2f}:{pts}", (">" if m < 0 else "<")
    return None


def default_heuristic(ctx):
    """Terminal heuristic — always fires; mirrors the historical final fallback."""
    return "default", "="


DEFAULT_HEURISTICS = [
    atp_synthase_heuristic,
    abc_transporter_heuristic,
    stored_bounds_heuristic,
    mmdeltag_band_heuristic,
    low_energy_heuristic,
    default_heuristic,
]


def make_ln_reversibility_index_heuristic(ln_ri_by_rxn, threshold=LN_RI_THRESHOLD):
    """Optional heuristic (not in ``DEFAULT_HEURISTICS``): eQuilibrator's
    ln-reversibility-index. ``ln_ri_by_rxn`` is ``{rxn_id: ln(gamma)}``.
    Serves as the reference for authoring energy-derived custom heuristics."""
    def heuristic(ctx):
        ln_ri = ln_ri_by_rxn.get(ctx.rxn_entry['id'])
        if ln_ri is not None and abs(ln_ri) > threshold:
            return f"lnRI: {ln_ri:.2f}", (">" if ln_ri < 0 else "<")
        return None
    return heuristic


# --- Pluggable energy sources: (rxn_entry) -> (dg, dge, source_label) -----
def top_level_energy(db_level):
    """Historical source: top-level ``deltag`` gated by ``db_level`` eligibility.
    ``db_level`` is ``''`` / ``'GC'`` / ``'EQ'`` / ``'DGP'``."""
    def resolve(rxn_entry):
        return _energy_for(rxn_entry, db_level)
    return resolve


def per_source_energy(label):
    """Source reading ``thermodynamics[label]``'s OWN dg (not top-level ``deltag``)."""
    def resolve(rxn_entry):
        pair = _thermo_pair(rxn_entry, label)
        if pair is None:
            return None, None, label
        return pair[0], pair[1], label
    return resolve


def explicit_energy(dg, dge):
    """Source wrapping an explicit ``(dg, dge)`` pair — for the per-source
    updaters and :func:`reversibility_from_energy`."""
    def resolve(rxn_entry):
        return dg, dge, None
    return resolve


# --- Cascade runner -------------------------------------------------------
def run_reversibility(rxn_entry, energy_source, heuristics=DEFAULT_HEURISTICS):
    """Resolve energy, then run ``heuristics`` until one fires. Returns
    ``(status, operator, source_label)``, or ``(None, None, label)`` when the
    energy source yields no usable energy (caller handles EMPTY/incomplete)."""
    dg, dge, source_label = energy_source(rxn_entry)
    if dg is None:
        return None, None, source_label
    ctx = Context(rxn_entry, float(dg), float(dge))
    for heuristic in heuristics:
        result = heuristic(ctx)
        if result is not None:
            return result[0], result[1], source_label
    return "default", "=", source_label
