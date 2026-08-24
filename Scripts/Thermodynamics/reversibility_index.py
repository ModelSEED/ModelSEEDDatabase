"""The eQuilibrator 2012 reversibility index, as a standalone callable.

Reference implementation of the reversibility index Gamma defined in

    Noor E, Bar-Even A, Flamholz A, Lubling Y, Davidi D, Milo R (2012)
    "An integrated open framework for thermodynamics of reactions that
    combines accuracy and coverage." Bioinformatics 28(15):2037-2044

and exposed by eQuilibrator 2.0 (Flamholz et al. 2012, NAR 40:D770) as
``ln_reversibility_index``. eQuilibrator 3.0 (Beber et al. 2022) still ships the
same quantity; it defines no index of its own.

Gamma answers one question: **by what fold change would every reactant
concentration have to move, symmetrically, before the reaction runs backwards?**
Gamma = 1000 is the canonical cut -- Noor 2012's window of 3 uM to 3 mM around a
100 uM reference -- so |ln Gamma| > ln(1000) marks a reaction that physiology
cannot plausibly reverse.

Definition
----------
For reagent coefficients ``nu_i`` (negative for substrates), with water and
protons dropped from both sums::

    dG'm  = dG'o + R T (sum nu_i) ln c            c = 1 mM
    ln G  = (2 / sum |nu_i|) * dG'm / (R T)

The ``2 / sum|nu|`` factor is what makes the index comparable across reactions
of different molecularity: it converts the free energy into the log-scale
concentration ratio a *single* reactant would need to span, given that
substrates move one way and products the other.

Conventions this module deliberately keeps
------------------------------------------
* **Water and protons are excluded** from both coefficient sums. This is
  eQuilibrator's convention (``Reaction.items(protons=False, water=False)``):
  the solvent is at unit activity and pH is a fixed condition of the transform,
  not a free concentration. ``exclude_protons=False`` is available for callers
  who want to see what a proton-translocating reaction looks like when the
  protons are allowed to count.
* **Every (compound, compartment) pair is its own species.** A reaction that
  moves a metabolite across a membrane is treated as ordinary chemistry -- two
  distinct species, one consumed and one produced -- with no membrane term, no
  compartment collapse, and no structural shortcut. This is the same rule
  applied to every reaction, which is the point: it is what lets a transport
  reaction be compared with a cytosolic one on the same axis. (It is also the
  opposite of the compartment-collapsing bug in
  ``Retrieve_eQuilibrator_Reactions_Energies.py``, which keys the MetaNetX
  formula on compound id alone and so cancels the two sides against each other.)

Units are kcal/mol throughout, matching what ModelSEED stores. Pass
``rt=RT_KJ`` to work in kJ/mol.

Usage
-----
    from reversibility_index import ln_reversibility_index, direction_from_index

    ln_gamma = ln_reversibility_index(rxn["stoichiometry"], dg)
    op = direction_from_index(ln_gamma)          # '>', '<' or '='
"""
from math import exp, inf, isfinite, log

__all__ = [
    "RT_KCAL", "RT_KJ", "PHYSIOLOGICAL_CONC", "LN_GAMMA_THRESHOLD",
    "PROTON", "WATER",
    "coefficient_sums", "physiological_dg_prime",
    "ln_reversibility_index", "ln_reversibility_index_error",
    "reversibility_index", "direction_from_index", "index_report",
]

# --- Constants ------------------------------------------------------------
TEMPERATURE = 298.15                       # K
GAS_CONSTANT_KCAL = 0.0019858775           # kcal / mol / K
RT_KCAL = TEMPERATURE * GAS_CONSTANT_KCAL  # 0.5921 kcal/mol
RT_KJ = RT_KCAL * 4.184                    # 2.4777 kJ/mol

# eQuilibrator's physiological convention: every aqueous reactant at 1 mM.
PHYSIOLOGICAL_CONC = 1.0e-3

# Noor 2012's default cut: a 1000-fold concentration swing.
LN_GAMMA_THRESHOLD = log(1000.0)           # 6.907755...

PROTON = "cpd00067"
WATER = "cpd00001"


# --- Core -----------------------------------------------------------------
def coefficient_sums(stoichiometry, exclude_protons=True, exclude_water=True,
                     proton_id=PROTON, water_id=WATER):
    """Return ``(sum_nu, sum_abs_nu, n_species)`` over the counted reagents.

    ``stoichiometry`` is ModelSEED's list of reagent dicts, each with at least
    ``compound`` and ``coefficient``; ``compartment`` is read when present so
    that the same compound on two sides of a membrane stays two species.

    Coefficients for repeated (compound, compartment) pairs are summed first, so
    a species written twice on the same side counts once with the combined
    coefficient -- and a species that genuinely cancels drops out.
    """
    merged = {}
    for rgt in stoichiometry:
        cpd = rgt["compound"]
        if exclude_water and cpd == water_id:
            continue
        if exclude_protons and cpd == proton_id:
            continue
        key = (cpd, rgt.get("compartment", 0))
        merged[key] = merged.get(key, 0.0) + float(rgt["coefficient"])

    sum_nu = sum_abs = 0.0
    n_species = 0
    for coeff in merged.values():
        if coeff == 0.0:
            continue
        sum_nu += coeff
        sum_abs += abs(coeff)
        n_species += 1
    return sum_nu, sum_abs, n_species


def physiological_dg_prime(standard_dg_prime, sum_nu, conc=PHYSIOLOGICAL_CONC,
                           rt=RT_KCAL):
    """dG'm -- the standard dG' moved to ``conc`` (1 mM) for every reagent.

    ``dG'm = dG'o + R T (sum nu) ln c``. Equivalent to
    ``ComponentContribution.physiological_dg_prime``.
    """
    return standard_dg_prime + rt * sum_nu * log(conc)


def ln_reversibility_index(stoichiometry, standard_dg_prime,
                           conc=PHYSIOLOGICAL_CONC, rt=RT_KCAL,
                           exclude_protons=True, exclude_water=True):
    """ln(Gamma) for one reaction. ``None`` when no reagents remain to count.

    ``standard_dg_prime`` is dG'o in the same energy unit as ``rt``
    (kcal/mol by default).
    """
    sum_nu, sum_abs, _ = coefficient_sums(
        stoichiometry, exclude_protons=exclude_protons,
        exclude_water=exclude_water)
    if sum_abs == 0.0:
        return None
    dgm = physiological_dg_prime(standard_dg_prime, sum_nu, conc, rt)
    return (2.0 / sum_abs) * dgm / rt


def ln_reversibility_index_error(stoichiometry, standard_dg_prime_error,
                                 exclude_protons=True, exclude_water=True,
                                 rt=RT_KCAL):
    """Uncertainty on ln(Gamma) propagated from the dG'o uncertainty.

    The concentration term is exact, so the error scales by the same
    ``2 / sum|nu| / RT`` factor as the energy itself.
    """
    _, sum_abs, _ = coefficient_sums(
        stoichiometry, exclude_protons=exclude_protons,
        exclude_water=exclude_water)
    if sum_abs == 0.0:
        return None
    return (2.0 / sum_abs) * abs(standard_dg_prime_error) / rt


def reversibility_index(stoichiometry, standard_dg_prime, **kwargs):
    """Gamma itself. Returns ``inf``/``0.0`` rather than overflowing, since
    ln(Gamma) routinely exceeds the float exponent range."""
    ln_gamma = ln_reversibility_index(stoichiometry, standard_dg_prime, **kwargs)
    if ln_gamma is None:
        return None
    if ln_gamma > 709.0:
        return inf
    if ln_gamma < -745.0:
        return 0.0
    return exp(ln_gamma)


def direction_from_index(ln_gamma, ln_gamma_err=0.0, z=0.0,
                         threshold=LN_GAMMA_THRESHOLD):
    """Map ln(Gamma) onto a ModelSEED direction operator.

    ``ln Gamma < -threshold`` means the forward direction is the spontaneous one
    (dG'm is negative), so the operator is ``>``; the mirror case gives ``<``;
    anything inside the window is ``=``.

    ``z`` is a confidence margin in units of the propagated sigma. ``z=0`` is the
    eQuilibrator 2.0 point-estimate rule as published; ``z=1`` additionally
    requires the one-sigma interval to clear the threshold.
    """
    if ln_gamma is None or not isfinite(ln_gamma):
        return "?"
    margin = abs(ln_gamma) - z * abs(ln_gamma_err or 0.0)
    if margin <= threshold:
        return "="
    return ">" if ln_gamma < 0 else "<"


def index_report(stoichiometry, standard_dg_prime, standard_dg_prime_error=0.0,
                 z=0.0, threshold=LN_GAMMA_THRESHOLD, **kwargs):
    """Everything the index says about one reaction, as a dict. Convenience
    wrapper for reports and tests."""
    sum_nu, sum_abs, n_species = coefficient_sums(
        stoichiometry,
        exclude_protons=kwargs.get("exclude_protons", True),
        exclude_water=kwargs.get("exclude_water", True))
    ln_gamma = ln_reversibility_index(stoichiometry, standard_dg_prime, **kwargs)
    err = ln_reversibility_index_error(
        stoichiometry, standard_dg_prime_error,
        exclude_protons=kwargs.get("exclude_protons", True),
        exclude_water=kwargs.get("exclude_water", True),
        rt=kwargs.get("rt", RT_KCAL))
    return {
        "sum_nu": sum_nu,
        "sum_abs_nu": sum_abs,
        "n_species": n_species,
        "dg_prime_m": (physiological_dg_prime(
            standard_dg_prime, sum_nu,
            kwargs.get("conc", PHYSIOLOGICAL_CONC),
            kwargs.get("rt", RT_KCAL)) if sum_abs else None),
        "ln_gamma": ln_gamma,
        "ln_gamma_err": err,
        "gamma": reversibility_index(stoichiometry, standard_dg_prime, **kwargs),
        "operator": direction_from_index(ln_gamma, err, z, threshold),
    }
