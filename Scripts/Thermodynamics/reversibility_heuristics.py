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

The gate that makes the EQ set differ most from GC exists because ``sigma`` of
1e5 kJ/mol is eQuilibrator's "could not decompose this reaction" marker, not an
error bar -- 4,934 reaction records carry it, and the GC bounds rule silently
swallows them (observed real sigma tops out at 65.35 kcal/mol, so the gap is
unambiguous).

TRANSPORT. There was a second gate, ``eq_transport_uncorrected_heuristic``,
returning "?" for every ``is_transport`` reaction. Removed 2026-09-08. It cited
two reasons and only one had survived:

  * The compartment-collapse defect is FIXED. It described the superseded
    MetaNetX-mediated retrieval, which keyed the formula on MetaNetX id and so
    discarded compartment, leaving 1,102 transport reactions carrying a dG for a
    different reaction. generate_modelseed_energies.py now sums coefficients
    across compartments deliberately -- a pure translocation cancels to nothing
    and is refused outright (4,078 such reactions carry no eQ value), while a
    driven pump keeps the chemistry that powers it.
  * The electrochemical term is still missing: Beber 2022 needs
    -N_H*RT*ln(10^dpH) - Q*F*dPhi across a membrane and this pipeline does not
    compute it. That is now a stated CAVEAT ON THE DATABASE rather than a
    refusal, matching how dGPredictor has always treated transport: every
    transport reaction is scored as ordinary biochemistry, without the membrane
    potential, and the reader is told so.

The rule was also barely reaching its cases -- of 1,687 eQ-scored transport
reactions, abc_transporter_heuristic decided 1,454 and atp_synthase_heuristic 15,
leaving it 218. Removing it moves 212 reactions out of "?": 121 to ">", 83 to
"=", 8 to "<".
"""
from math import log

# Single implementation of the Noor 2012 / eQuilibrator 2.0 reversibility index.
# ``Context.ln_gamma`` delegates to it so the cascade and any standalone caller
# can never drift apart.
from reversibility_index import (
    LN_GAMMA_THRESHOLD, PHYSIOLOGICAL_CONC,
    coefficient_sums, direction_from_index,
    ln_reversibility_index, ln_reversibility_index_error,
)

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

# A predecessor dGPredictor record, broken by RDKit canonicalization drift, was
# retired in Aug 2026 and replaced by the retrained model that ships as
# "dGPredictor". Its name and its DGPM report level were removed 2026-09-08;
# nothing in the shipped Biochemistry/*.json ever carried them.
DB_LEVEL_LABEL = {"GC": "Group contribution", "EQ": "eQuilibrator",
                  "DGP": "dGPredictor"}
DB_LEVEL_NOTE = {"GC": "GCC", "EQ": "EQU"}
DB_LEVEL_PRIORITY = ("EQ", "GC", "DGP")

# Levels whose energy must be read from ``thermodynamics[label]`` rather than the
# flat top-level ``deltag``/``deltagerr``. The flat field is written by whichever
# source last promoted to canonical, so for every source but that one it holds a
# *different* source's number under this source's name. GC is deliberately absent:
# its historical report is byte-compared against upstream, and the flat field is
# in fact the GC value for the reactions that report covers.
PER_SOURCE_LEVELS = ("EQ", "DGP")

LN_RI_THRESHOLD = LN_GAMMA_THRESHOLD    # ln(1000); Noor 2012 default

# --- eQuilibrator-specific constants --------------------------------------
KJ_PER_KCAL = 4.184

# eQuilibrator reports sigma = 1e5 kJ/mol (= 23,900.57 kcal/mol) for a reaction
# it cannot decompose; multiples of it appear when several degrees of freedom are
# unknown. This is a MARKER, not an error bar.
#
# THE CUT IS NOT THE MARKER. It sits an order of magnitude below it, in the empty
# gap that separates real uncertainties from refusals, so that a marker scaled
# down by propagation is still caught. Measured over the shipped
# Biochemistry/reaction_*.json on 2026-09-08:
#
#   largest GENUINE eQuilibrator sigma        40.24 kcal/mol
#   ---- the gap ----
#   smallest MARKER value observed        12,482.48 kcal/mol  (= 52,227 kJ/mol)
#
# The smallest observed marker is about half the nominal 1e5 kJ/mol, which is why
# the cut is not placed at the nominal value. An earlier version of this comment
# quoted 65.35 and 23,900.57 for those bounds; those came from
# MetaNetX_Reaction_Energies.tbl, a different corpus from the shipped reactions
# the gate actually runs against.
#
# Set to a round 2,500 kcal/mol on request 2026-09-08, replacing 1e4/KJ_PER_KCAL
# (2390.06). Both sit in the same gap and the change decides 0 reactions; 2,500
# is 62x the largest genuine sigma and 5.0x below the smallest marker. It also
# clears dGPredictor's largest sigma (1,427.33) by 1.75x, which matters now that
# the gate is shared with that source.
EQ_UNDECOMPOSABLE_SIGMA = 2500.0                     # kcal/mol

# eQuilibrator's physiological convention: every aqueous reactant at 1 mM.
# dG'm = dG'^o + RT * sum(nu) * ln(EQ_PHYSIOLOGICAL_CONC), water and protons
# excluded from the sum (equilibrator_cache.reaction.items(protons=False,
# water=False)). Note this is flat, unlike the GC cascade's per-compound
# CELL_CONC / CO2 / O2 / H2 concentrations.
EQ_PHYSIOLOGICAL_CONC = PHYSIOLOGICAL_CONC

# Confidence margin, in units of the propagated ln(Gamma) sigma, that the
# reversibility index must clear before the reaction is called irreversible.
# 1.0 for eQuilibrator 3.0 (Beber 2022 makes uncertainty first-class); 0.0
# reproduces the eQuilibrator 2.0 point-estimate behaviour.
EQ_CONFIDENCE_Z = 1.0

# Fraction of the ln(1000) threshold below which the propagated ln(Gamma) sigma
# is too small for the confidence margin to mean anything. 0.15 x ln(1000) =
# 1.04. Below that the error bar cannot move the reaction across the threshold
# in any meaningful sense, so the point estimate decides and "?" is not used.
RI_TIGHT_FRACTION = 0.15


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
        sum_nu, _, _ = coefficient_sums(self.rxn_entry['stoichiometry'])
        return self.dg + RT_CONST * sum_nu * log(EQ_PHYSIOLOGICAL_CONC)

    @property
    def ln_gamma(self):
        """Reversibility index in natural log, ``(2 / sum|nu|) * dG'm / RT``
        (Noor et al. 2012). ``None`` when the reaction has no reagents left
        after dropping water and protons, mirroring eQuilibrator's guard.

        Interpretation: ``Gamma`` is the fold change every reactant
        concentration must undergo to reverse the reaction, so the sign follows
        dG'm and the magnitude says how hard the reversal is.

        Delegates to :func:`reversibility_index.ln_reversibility_index` -- see
        that module for the conventions, including why every
        (compound, compartment) pair counts as its own species."""
        if self._ln_gamma is None:
            self._ln_gamma = ln_reversibility_index(
                self.rxn_entry['stoichiometry'], self.dg, rt=RT_CONST)
        return self._ln_gamma

    @property
    def ln_gamma_err(self):
        """``ln_gamma`` propagated from the reported dG uncertainty. The
        concentration term is exact, so only ``dge`` carries through."""
        return ln_reversibility_index_error(
            self.rxn_entry['stoichiometry'], self.dge, rt=RT_CONST)


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


def make_sentinel_heuristic(sigma_gate=None):
    """Reject a source's "no estimate" marker before any rule reads it.

    ONE rule for one concept: the source published a number and disowned it in
    the same breath. Sources differ only in which field carries the message --
    Group Contribution puts it in the energy (dg = 1e7), eQuilibrator in the
    uncertainty (a cut at 2,500 kcal/mol, below its 1e5 kJ/mol marker)
    -- so both are checked here
    rather than in two rules at different depths, which is what they were until
    2026-09-08.

    ``sigma_gate`` is the "predictor declined" cut, or None for a source with no
    such marker. dGPredictor passes None deliberately: a fragment the model has
    never seen contributes nothing to the sum, so there is no sentinel to gate
    on and silent extrapolation is the failure mode instead of a loud refusal.

    Placing this first in every rule set makes the guarantee hold regardless of
    entry point. The dg check previously lived only OUTSIDE the cascade, in
    _thermo_pair and reversibility_from_energy, so a caller reaching
    run_reversibility through explicit_energy skipped it and the sentinel flowed
    into stored_bounds_heuristic, which returned a confident "<" off
    "MdeltaG(Min): 9999994.50".
    """
    def sentinel_energy_heuristic(ctx, _cut=sigma_gate):
        if ctx.dg == SENTINEL_DG:
            return f"no estimate: dg sentinel {SENTINEL_DG:.0f}", "?"
        if _cut is not None and abs(ctx.dge) >= _cut:
            return f"no estimate: sigma {ctx.dge:.0f}", "?"
        return None
    return sentinel_energy_heuristic


# Default instance for sources whose only marker is the dg sentinel.
sentinel_energy_heuristic = make_sentinel_heuristic()


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
    """Terminal rule for GC: NOTHING ABOVE COULD DECIDE, so report "?".

    Returned "=" until 2026-09-08, which asserted reversibility for 10,201
    reactions on the strength of no test at all -- their median sigma is 12.24
    kcal/mol and 70% of them have sigma >= |mMdeltaG|, so the energy's own SIGN
    is not established, let alone its magnitude. "?" is what the cascade
    actually knows here.
    """
    return "default", "?"


def canonical_default_heuristic(ctx):
    """Terminal rule for the CANONICAL top-level ``reversibility`` field: "=".

    THE SCOPE SPLIT, 2026-09-08. default_heuristic was changed to "?" so the
    per-source thermodynamics dict stops asserting reversibility it cannot
    support. That change was scoped to the thermodynamics field ONLY. It leaked
    into the canonical field because both paths shared one GC rule list, and
    Apply_2020_Reversibility_Policy.py runs the GC set -- so regenerating the
    canonical field silently moved 8,687 reactions from "=" to "?".

    The canonical field is the 2020 series and must stay comparable across
    releases, exactly as GC's ABC-transporter rule is kept for continuity. It
    keeps the historical "=" terminal here while GC_HEURISTICS keeps "?" for
    per-source use. Selected by ENTRY POINT, not caller discipline:
    get_heuristics() serves the canonical path, heuristics_for_source() the
    per-source path.
    """
    return "default", "="


GC_HEURISTICS = [
    sentinel_energy_heuristic,
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

    Prefer the ``ri_index_heuristic`` built by :func:`make_ri_heuristics`,
    which derives ln(Gamma)
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

def make_ri_heuristics(z=0.0, threshold=LN_RI_THRESHOLD,
                       sigma_gate=EQ_UNDECOMPOSABLE_SIGMA):
    """The Noor 2012 reversibility index, for any source publishing dG + sigma.

    THE ONE CASCADE FOR BOTH eQuilibrator AND dGPredictor. Until 2026-09-08
    there were two, ``make_eq_heuristics`` and this one, whose index and default
    rules were byte-for-byte equivalent apart from an ``EQ:``/``RI:`` label
    prefix and a defensive ``or 0.0`` on the sigma. Measured across all 110,703
    scorable contexts they returned the SAME operator every time, and
    ``ln_gamma_err`` was never None, so the guard never fired. Two names for one
    heuristic is how they drift apart, so the EQ copies were deleted and both
    sources now build from here. Status strings are prefixed ``RI:`` for both;
    the source is carried separately by ``run_reversibility``, so nothing is lost.

    Carries exactly one structural rule, ATP synthase -- see below. It carries
    NO ABC-transporter shortcut and no membrane gate: a reaction that moves a
    species across a compartment is scored like any other, from its own
    stoichiometry and its own dG. That puts every reaction on one axis, so a
    transport call can be compared with a cytosolic one instead of being decided
    by a rule that fires before the energy is ever read.

    ``sigma_gate`` is the "predictor declined" cut. Pass ``None`` to disable it
    for a source that has no such marker.
    """
    # Sentinel first, for the same reason as GC: a caller reaching the cascade
    # through explicit_energy skips the outer guards entirely. One rule,
    # carrying this source's sigma marker if it has one.
    #
    # Then ATP SYNTHASE, THE ONE STRUCTURAL RULE THIS CASCADE NEEDS. Without it
    # the index calls these confidently and wrongly: lnGamma = +/-11.70 +/- 0.09,
    # 1.7x the ln(1000) threshold, so no uncertainty gate can soften it. The 15
    # records are the same chemistry written both ways -- 7 as hydrolysis, 8 as
    # synthesis -- so the index returns 7 ">" and 8 "<", tracking transcription
    # order rather than biology. The proton-motive force is the entire driving
    # force here and this pipeline does not compute it, so "=" is the only
    # defensible call. Applied to dGPredictor too from 2026-09-08: the exposure
    # is identical and was simply never put to that source.
    #
    # NOT abc_transporter_heuristic. Measured against the 1,454 reactions it had
    # decided under EQ: the index reaches the same answer for 1,414 with ZERO
    # reversals, so the ATP hydrolysis in the stoichiometry already drives the
    # call. Of the 40 differences, 26 contain no ATP hydrolysis at all -- the
    # rule keys on phosphate-count sign and had misidentified them -- and 14 are
    # genuine ATP-driven transport sitting inside the reversible band, now "=".
    # GC keeps it: that cascade is Chris's historical one, preserved as history.
    rules = [make_sentinel_heuristic(sigma_gate), atp_synthase_heuristic]

    def index_rule(ctx, _z=z, _thr=threshold):
        ln_gamma = ctx.ln_gamma
        if ln_gamma is None:
            return None
        err = ctx.ln_gamma_err or 0.0
        if abs(ln_gamma) - _z * err > _thr:
            return (f"RI:lnGamma: {ln_gamma:.2f}+/-{err:.2f}",
                    ">" if ln_gamma < 0 else "<")
        return None
    index_rule.__name__ = 'ri_index_heuristic'
    rules.append(index_rule)

    def terminal(ctx, _z=z, _thr=threshold):
        ln_gamma = ctx.ln_gamma
        if ln_gamma is None:
            return "RI:no-reagents", "="
        err = ctx.ln_gamma_err or 0.0
        g = abs(ln_gamma)
        # "=" MEANS REVERSIBLE and "?" MEANS UNKNOWN. Three states reach here:
        #
        #  1. The error bar fits INSIDE the band -> positively reversible, "=".
        #  2. The bar STRADDLES ln(1000) but is small relative to it: the
        #     reaction is sitting ON the threshold with a well-determined
        #     energy. Calling this "?" said "no evidence" about reactions whose
        #     energy is among the best in the database (rxn01101: dG = 8.09
        #     +/- 0.14, |lnGamma| = 6.83 +/- 0.12 against a threshold of 6.91).
        #     It reports "=", NOT a direction, and deliberately so. Straddling
        #     means |g - thr| <= err BY DEFINITION, so the point estimate is
        #     always within one sigma of the line -- shrinking sigma moves the
        #     line closer in absolute terms rather than resolving which side the
        #     reaction is on. A hard ">" here would be a ~65%-confidence call
        #     dressed as a determination. "=" is the permissive reading that the
        #     evidence does support: at physiological concentrations this
        #     reaction is marginal, and concentration control can cross it.
        #  3. The bar straddles and is wide -> genuinely unknown, "?".
        #
        # Rejected: applying the point estimate unconditionally. On dGPredictor,
        # whose median propagated sigma is 15.51 against a threshold of 6.91,
        # that would hand out 11,251 hard directional calls read off noise.
        # The gate is what makes the point estimate defensible.
        #
        # Until 2026-09-08 case 2 shipped as "?" and, before that, cases 2 and 3
        # both shipped as "=".
        if g + _z * err < _thr:
            return f"RI:reversible: {ln_gamma:.2f}+/-{err:.2f}", "="
        if err < RI_TIGHT_FRACTION * _thr:
            return f"RI:near-threshold: {ln_gamma:.2f}+/-{err:.2f}", "="
        return f"RI:unknown: {ln_gamma:.2f}+/-{err:.2f}", "?"
    terminal.__name__ = 'ri_default_heuristic'
    rules.append(terminal)
    return rules


# eQuilibrator and both dGPredictor sources now build from ONE factory with
# IDENTICAL arguments. Keeping two constants preserves the registry seam (and the
# EQ/DGP report levels) should they ever need to diverge again; today they
# do not, and any change must be made to both deliberately rather than to one by
# accident.
#
# THE SIGMA GATE IS NOW SHARED, and the separation is clean. dGPredictor has no
# "could not decompose" marker -- an unseen fragment contributes nothing to the
# sum, so silent extrapolation is its failure mode and no sigma cut catches it --
# but the cut is carried anyway, for consistency and as a tripwire. Measured
# 2026-09-08 over every shipped thermodynamics entry:
#
#   source                  n       median      max sigma   >= 2390.06
#   dGPredictor          29,617      17.01        1,427.33        0
#   eQuilibrator         25,175       0.77      861,955.72    3,386   <- markers
#   Group contribution   56,002      24.61   10,000,000.00   26,555   <- markers
#
# Excluding each source's marker records, the largest GENUINE sigma anywhere is
# 566.61 (GC), and eQuilibrator's real values stop at 40.24 before jumping to the
# refusal scale. So 2,390.06 sits above every real uncertainty in the database
# with room to spare, and applying it to dGPredictor changes 0 reactions today.
# Do NOT read it as evidence that dGPredictor declines -- it does not.
#
# GC is deliberately NOT given the gate: it encodes its marker in BOTH fields, so
# 26,555 of its entries would trip a sigma cut, and the dg == 1e7 test in the
# same rule already catches every one of them.
#
# THE ONE-SIGMA MARGIN COSTS THE TWO SOURCES VERY DIFFERENTLY, and that is the
# honest reading of their error bars rather than a defect. Threshold = ln(1000)
# = 6.91. Measured 2026-09-08 on the propagated ln(Gamma) sigma, which is what
# the index actually consumes:
#
#   source          n      med |lnG|   med sigma_lnG   z=0 directional -> z=1
#   eQuilibrator  21,789       6.77            0.62      10,814 -> 10,209  (-5.6%)
#   dGPredictor   29,616      11.60           15.51      20,145 ->  8,894 (-55.9%)
#
# eQuilibrator's typical error bar is 11x BELOW the decision threshold;
# dGPredictor's is 2.2x ABOVE it. So dGPredictor makes bolder point predictions
# (median |lnGamma| 11.60 vs 6.77) and has far less right to them: taken at face
# value it would out-call eQuilibrator nearly two to one, and once the error bars
# are honoured it falls below it. The margin is kept precisely because silent
# extrapolation is its failure mode and nothing else catches it.
#
# An earlier version of this comment claimed dGPredictor "reports a fit residual
# (median 0.35, max 6.1 kcal/mol)" and that z=1 moved only "784 of 27,715"
# reactions to ambiguous. Both are wrong by more than an order of magnitude: the
# shipped raw sigma has median 17.01 and max 1,427.33, and z=1 sends 18,657 of
# 29,616 to "?". Do not restore those figures.
#
# The point-estimate answer is available by building the set directly with
# make_ri_heuristics(z=0.0); it is no longer registered. The RI set was removed
# 2026-09-08 alongside EQ2 -- both were comparison arms production never
# selected, and a registry entry is an invitation to select one.
#
EQ_HEURISTICS = make_ri_heuristics(z=EQ_CONFIDENCE_Z,
                                   sigma_gate=EQ_UNDECOMPOSABLE_SIGMA)
DGP_HEURISTICS = make_ri_heuristics(z=EQ_CONFIDENCE_Z,
                                    sigma_gate=EQ_UNDECOMPOSABLE_SIGMA)


# --- Rule-set registry ----------------------------------------------------
# GC_CANONICAL is GC with the historical "=" terminal, for the top-level field.
GC_CANONICAL_HEURISTICS = GC_HEURISTICS[:-1] + [canonical_default_heuristic]

HEURISTIC_SETS = {
    'GC': GC_HEURISTICS,
    'GC_CANONICAL': GC_CANONICAL_HEURISTICS,
    'EQ': EQ_HEURISTICS,
    'DGP': DGP_HEURISTICS,
}

DEFAULT_HEURISTIC_SET = 'GC'

# ``thermodynamics`` subkey -> rule-set name. Only genuinely UNKNOWN labels fall
# back to GC; every source shipped today is mapped explicitly below. (The old
# comment here claimed the dGPredictor sources wanted the GC fallback -- stale
# since DGP was added, and contradicted by the dict immediately below it.)
SOURCE_HEURISTIC_SET = {
    'eQuilibrator': 'EQ',
    'dGPredictor': 'DGP',
}


def get_heuristics(name=None):
    """Rule list for the CANONICAL top-level path. Unknown/missing name -> GC.

    Returns the GC_CANONICAL variant wherever plain GC is asked for, because
    every caller of this function writes the top-level ``reversibility`` field
    and that field keeps the historical "=" terminal. The per-source path uses
    heuristics_for_source() instead. See canonical_default_heuristic.
    """
    resolved = HEURISTIC_SETS.get(name or DEFAULT_HEURISTIC_SET, GC_HEURISTICS)
    return GC_CANONICAL_HEURISTICS if resolved is GC_HEURISTICS else resolved


def heuristic_set_for_source(label=None):
    """Rule-set *name* appropriate to a ``thermodynamics`` subkey."""
    return SOURCE_HEURISTIC_SET.get(label, DEFAULT_HEURISTIC_SET)


def heuristics_for_source(label=None):
    """Rule *list* for a ``thermodynamics`` subkey -- the PER-SOURCE path.

    Reads HEURISTIC_SETS directly rather than via get_heuristics(), which
    redirects plain GC to the canonical "=" terminal for the top-level field.
    This path wants the "?" terminal. GC by default.
    """
    name = heuristic_set_for_source(label)
    return HEURISTIC_SETS.get(name or DEFAULT_HEURISTIC_SET, GC_HEURISTICS)


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

    ``DGP`` reads its own sublist for the same reason: since the
    additive-thermodynamics refactor nothing overwrites ``deltag``, so scoring a
    dGPredictor level off the flat field scores the Group-Contribution number
    and labels it dGPredictor.

    ``GC`` and the unfiltered run keep the historical top-level source."""
    if db_level in PER_SOURCE_LEVELS:
        return per_source_energy(DB_LEVEL_LABEL[db_level])
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
    # Unreachable: every registered set ends in a terminal rule. Kept as a
    # backstop, and "?" for the same reason default_heuristic returns it -- a
    # cascade that decided nothing has not established reversibility.
    return "default", "?", source_label
