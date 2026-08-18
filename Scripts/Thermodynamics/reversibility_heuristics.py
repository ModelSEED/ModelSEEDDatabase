"""Composable reaction-reversibility cascade: pluggable heuristics + energy sources.

The cascade is (energy_source, [heuristics]). The energy source produces
``(dg, dge, source_label)``; each heuristic is a ``(ctx: Context) -> (status, op)
| None`` callable, and the first non-None result wins.

Extend by:
  - **New energy source**: write ``(rxn_entry) -> (dg, dge, label)``.
    See :func:`top_level_energy` / :func:`per_source_energy` / :func:`explicit_energy`.
  - **New heuristic**: write ``(ctx) -> (status, op) | None`` and insert into
    the list passed to :func:`run_reversibility`, e.g.
    ``GC_HEURISTICS + [make_ln_reversibility_index_heuristic(ln_ri)]``.
  - **New rule set**: add it to :data:`HEURISTIC_SETS` (and, if it belongs to a
    particular ``thermodynamics`` subkey, to :data:`SOURCE_HEURISTIC_SET`).

Rule sets are selected **per thermodynamic data source**, because the sources
do not fail in the same way:

``GC_HEURISTICS`` (default)
    The historical Group-Contribution cascade, reproduced byte-for-byte (same
    order, same status strings). The regression test byte-compares the generated
    reversibility report against upstream/dev, so any change here needs to
    preserve those exact status strings and operator returns.

``EQ_HEURISTICS`` (eQuilibrator 3.0; Beber et al. 2022, NAR 50:D603)
    Built on eQuilibrator's own methodology rather than the GC concentration
    bounds. Directionality comes from the reversibility index of Noor et al.
    2012 (Bioinformatics 28:2037) as used by eQuilibrator 2.0 (Flamholz et al.
    2012, NAR 40:D770) -- Beber 2022 itself defines no reversibility index. What
    3.0 contributes, and what this set adds on top of the 2.0 rule, is a serious
    treatment of uncertainty: eQuilibrator's sigma gates the call, and reactions
    whose estimate it cannot stand behind return "?" rather than a permissive
    "=".

``EQ2_HEURISTICS`` (eQuilibrator 2.0; Flamholz 2012 + Noor 2012)
    The same cascade with the confidence margin switched off, i.e. the ln(Gamma)
    point estimate alone. Kept for comparison against the 3.0 behaviour.

The two gates that make the EQ set differ most from GC exist because the stored
eQuilibrator numbers have two known failure modes (see the README):
  * ``sigma`` of 1e5 kJ/mol is eQuilibrator's "could not decompose this reaction"
    marker, not an error bar -- 4,934 reaction records carry it, and the GC
    bounds rule silently swallows them (observed real sigma tops out at 65.35
    kcal/mol, so the gap is unambiguous).
  * ``Retrieve_eQuilibrator_Reactions_Energies.py`` keys its reaction formula on
    MetaNetX id and so discards compartment, collapsing any species present on
    both sides. 1,102 transport reactions therefore carry a dG for a *different*
    reaction. Beber 2022 additionally notes the transformed-ensemble framework
    is invalid across membranes without the -N_H*RT*ln(10^dpH) - Q*F*dPhi term,
    which this pipeline never applies.
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

# --- eQuilibrator-specific constants --------------------------------------
KJ_PER_KCAL = 4.184

# eQuilibrator reports sigma = 1e5 kJ/mol for a reaction it cannot decompose;
# multiples of it appear when several degrees of freedom are unknown. This is a
# marker, not an error bar. The cut sits in the empty gap between the largest
# genuine sigma observed in MetaNetX_Reaction_Energies.tbl (65.35 kcal/mol) and
# the smallest marker value (1e5 kJ/mol = 23900.57 kcal/mol).
EQ_UNDECOMPOSABLE_SIGMA = 1.0e4 / KJ_PER_KCAL        # 2390.06 kcal/mol

# eQuilibrator's physiological convention: every aqueous reactant at 1 mM.
# dG'm = dG'^o + RT * sum(nu) * ln(EQ_PHYSIOLOGICAL_CONC), water and protons
# excluded from the sum (equilibrator_cache.reaction.items(protons=False,
# water=False)). Note this is flat, unlike the GC cascade's per-compound
# CELL_CONC / CO2 / O2 / H2 concentrations.
EQ_PHYSIOLOGICAL_CONC = 0.001

# Confidence margin, in units of the propagated ln(Gamma) sigma, that the
# reversibility index must clear before the reaction is called irreversible.
# 1.0 for eQuilibrator 3.0 (Beber 2022 makes uncertainty first-class); 0.0
# reproduces the eQuilibrator 2.0 point-estimate behaviour.
EQ_CONFIDENCE_Z = 1.0


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
    nu_sum = abs_nu_sum = 0.0
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
        # eQuilibrator's two coefficient sums, over the same water/proton-free
        # reagent set: net (concentration correction) and absolute (ln Gamma).
        nu_sum += coeff
        abs_nu_sum += abs(coeff)
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
            'nu_sum': nu_sum, 'abs_nu_sum': abs_nu_sum,
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
    __slots__ = ('rxn_entry', 'dg', 'dge', '_terms', '_mMdeltaG', '_ln_gamma')

    def __init__(self, rxn_entry, dg, dge):
        self.rxn_entry, self.dg, self.dge = rxn_entry, dg, dge
        self._terms = None
        self._mMdeltaG = None
        self._ln_gamma = None

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

    @property
    def dg_prime_m(self):
        """eQuilibrator's physiological dG'm: every aqueous reactant at 1 mM.

        The GC cascade's :attr:`mMdeltaG` is the same idea with per-compound
        local concentrations (CO2 at 0.1 mM, O2/H2 at 1 uM); this one is flat,
        matching ``ComponentContribution.physiological_dg_prime``."""
        return self.dg + RT_CONST * self.terms['nu_sum'] * log(EQ_PHYSIOLOGICAL_CONC)

    @property
    def ln_gamma(self):
        """Reversibility index in natural log, ``(2 / sum|nu|) * dG'm / RT``
        (Noor et al. 2012). ``None`` when the reaction has no reagents left
        after dropping water and protons, mirroring eQuilibrator's guard.

        Interpretation: ``Gamma`` is the fold change every reactant
        concentration must undergo to reverse the reaction, so the sign follows
        dG'm and the magnitude says how hard the reversal is."""
        if self._ln_gamma is None:
            abs_nu = self.terms['abs_nu_sum']
            if abs_nu == 0:
                return None
            self._ln_gamma = (2.0 / abs_nu) * self.dg_prime_m / RT_CONST
        return self._ln_gamma

    @property
    def ln_gamma_err(self):
        """``ln_gamma`` propagated from the reported dG uncertainty. The
        concentration term is exact, so only ``dge`` carries through."""
        abs_nu = self.terms['abs_nu_sum']
        if abs_nu == 0:
            return None
        return (2.0 / abs_nu) * abs(self.dge) / RT_CONST


# --- Shared / Group-Contribution heuristics -------------------------------
# Signature: (ctx: Context) -> (status_label, operator) | None
# First non-None wins. Extend by defining another and appending to a rules list.
#
# ``atp_synthase_heuristic`` and ``abc_transporter_heuristic`` are structural,
# not energy-derived, so both the GC and EQ rule sets reuse them as-is.

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


GC_HEURISTICS = [
    atp_synthase_heuristic,
    abc_transporter_heuristic,
    stored_bounds_heuristic,
    mmdeltag_band_heuristic,
    low_energy_heuristic,
    default_heuristic,
]

# Back-compat alias: GC remains the default rule set for every source that has
# no set of its own. Existing importers keep working unchanged.
DEFAULT_HEURISTICS = GC_HEURISTICS


def make_ln_reversibility_index_heuristic(ln_ri_by_rxn, threshold=LN_RI_THRESHOLD):
    """Heuristic driven by a precomputed ``{rxn_id: ln(gamma)}`` map, e.g. the
    fourth column of ``eQuilibrator/MetaNetX_Reaction_Energies.tbl``.

    Prefer :func:`eq_reversibility_index_heuristic`, which derives ln(Gamma)
    from the stored dG and the reaction's own stoichiometry and so stays correct
    for the reactions where eQuilibrator scored a compartment-collapsed formula.
    Kept for callers that want to inject eQuilibrator's own published values."""
    def heuristic(ctx):
        ln_ri = ln_ri_by_rxn.get(ctx.rxn_entry['id'])
        if ln_ri is not None and abs(ln_ri) > threshold:
            return f"lnRI: {ln_ri:.2f}", (">" if ln_ri < 0 else "<")
        return None
    return heuristic


# --- eQuilibrator heuristics ----------------------------------------------
# Beber et al. 2022 (eQuilibrator 3.0) for the uncertainty treatment; Noor et
# al. 2012 / Flamholz et al. 2012 (eQuilibrator 2.0) for the reversibility
# index that supplies the actual direction.

def eq_undecomposable_heuristic(ctx):
    """eQuilibrator could not decompose the reaction, and says so with a
    ~1e5 kJ/mol sigma. There is no information in the accompanying dG, so
    report "?" rather than let a meaningless number reach the index rule."""
    if abs(ctx.dge) >= EQ_UNDECOMPOSABLE_SIGMA:
        return f"EQ:undecomposable: {ctx.dge:.0f}", "?"
    return None


def eq_transport_uncorrected_heuristic(ctx):
    """Transport reaction whose energy we cannot trust.

    Two independent reasons, both documented in the module docstring: the
    retrieval step collapses compartments when it builds the MetaNetX formula,
    and Beber 2022 notes the transformed framework needs a
    ``-N_H*RT*ln(10^dpH) - Q*F*dPhi`` term across a membrane that this pipeline
    never applies. ATP synthase and ABC transporters are decided structurally
    before this rule, so they never reach it."""
    if ctx.rxn_entry.get('is_transport') == 1:
        return "EQ:transport-uncorrected", "?"
    return None


def make_eq_reversibility_index_heuristic(z=EQ_CONFIDENCE_Z,
                                          threshold=LN_RI_THRESHOLD):
    """Direction from the reversibility index, requiring ``z`` sigma of margin.

    ``|ln Gamma| - z*sigma > ln(1000)`` means even the pessimistic end of the
    interval needs more than a 1000-fold concentration swing to reverse the
    reaction -- Noor 2012's headline window of 3 uM to 3 mM around 100 uM.
    ``z=0`` reduces this to the eQuilibrator 2.0 point-estimate test."""
    def heuristic(ctx):
        ln_gamma = ctx.ln_gamma
        if ln_gamma is None:
            return None
        margin = abs(ln_gamma) - z * ctx.ln_gamma_err
        if margin > threshold:
            return (f"EQ:lnGamma: {ln_gamma:.2f}+/-{ctx.ln_gamma_err:.2f}",
                    ">" if ln_gamma < 0 else "<")
        return None
    heuristic.__name__ = 'eq_reversibility_index_heuristic'
    return heuristic


def make_eq_default_heuristic(z=EQ_CONFIDENCE_Z, threshold=LN_RI_THRESHOLD):
    """Terminal EQ rule -- always fires, always "=".

    Splits the label so the report distinguishes a reaction that is confidently
    inside the reversible window from one whose interval merely straddles the
    threshold. Both are called reversible: Noor 2012 treats Gamma as a
    continuous index and reserves the directional call for clear cases, and "="
    is what the GC cascade's terminal rule returns too."""
    def heuristic(ctx):
        ln_gamma = ctx.ln_gamma
        if ln_gamma is None:
            return "EQ:no-reagents", "="
        err = ctx.ln_gamma_err
        state = "reversible" if abs(ln_gamma) + z * err < threshold else "ambiguous"
        return f"EQ:{state}: {ln_gamma:.2f}+/-{err:.2f}", "="
    heuristic.__name__ = 'eq_default_heuristic'
    return heuristic


def make_eq_heuristics(z=EQ_CONFIDENCE_Z, threshold=LN_RI_THRESHOLD):
    """Assemble an eQuilibrator rule set. ``z`` is the confidence margin in
    units of the propagated ln(Gamma) sigma (1.0 = eQuilibrator 3.0, 0.0 =
    eQuilibrator 2.0)."""
    return [
        # Structural first: these need no energy, so the two eQuilibrator data
        # defects cannot reach them.
        atp_synthase_heuristic,
        abc_transporter_heuristic,
        eq_undecomposable_heuristic,
        eq_transport_uncorrected_heuristic,
        make_eq_reversibility_index_heuristic(z, threshold),
        make_eq_default_heuristic(z, threshold),
    ]


# Module-level singletons so ``is`` comparisons in tests and callers are stable.
eq_reversibility_index_heuristic = make_eq_reversibility_index_heuristic()
eq_default_heuristic = make_eq_default_heuristic()

EQ_HEURISTICS = make_eq_heuristics()             # eQuilibrator 3.0, Beber 2022
EQ2_HEURISTICS = make_eq_heuristics(z=0.0)       # eQuilibrator 2.0, Flamholz 2012


# --- Rule-set registry ----------------------------------------------------
HEURISTIC_SETS = {
    'GC': GC_HEURISTICS,
    'EQ': EQ_HEURISTICS,
    'EQ2': EQ2_HEURISTICS,
}

DEFAULT_HEURISTIC_SET = 'GC'

# ``thermodynamics`` subkey -> rule-set name. Anything absent falls back to GC,
# which is what the dGPredictor sources want: they are bare dG predictions with
# no eQuilibrator-style uncertainty semantics behind them.
SOURCE_HEURISTIC_SET = {
    'eQuilibrator': 'EQ',
}


def get_heuristics(name=None):
    """Rule list for a rule-set name. Unknown or missing name -> GC."""
    return HEURISTIC_SETS.get(name or DEFAULT_HEURISTIC_SET, GC_HEURISTICS)


def heuristic_set_for_source(label=None):
    """Rule-set *name* appropriate to a ``thermodynamics`` subkey."""
    return SOURCE_HEURISTIC_SET.get(label, DEFAULT_HEURISTIC_SET)


def heuristics_for_source(label=None):
    """Rule *list* appropriate to a ``thermodynamics`` subkey. GC by default."""
    return get_heuristics(heuristic_set_for_source(label))


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


def energy_source_for_level(db_level):
    """Energy source that pairs naturally with a ``db_level``.

    ``EQ`` reads the eQuilibrator sublist's own dG and sigma rather than the
    top-level ``deltag``. That top-level value is only the eQuilibrator estimate
    for 1,797 of the 25,028 reactions that have one -- since the additive-
    thermodynamics refactor no caller overwrites ``deltag``, so the EQ run was
    scoring the Group-Contribution number and labelling it eQuilibrator. The EQ
    rule set also needs eQuilibrator's own sigma for its undecomposable gate,
    which the top-level ``deltagerr`` never carries.

    Every other level keeps the historical top-level source."""
    if db_level == 'EQ':
        return per_source_energy(DB_LEVEL_LABEL['EQ'])
    return top_level_energy(db_level)


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
