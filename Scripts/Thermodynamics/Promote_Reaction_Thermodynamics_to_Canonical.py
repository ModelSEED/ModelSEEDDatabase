#!/usr/bin/env python
"""Promote stranded per-source thermodynamics into the canonical deltag.

After the "additive descriptions" refactor, the Update_Reaction_*_Energies.py
scripts write only into each reaction's additive `thermodynamics` dict; none
write the canonical top-level deltag/deltagerr anymore. As a result many
reactions carry a perfectly good computed energy in `thermodynamics` while their
canonical deltag is still the 10000000 sentinel (and reversibility is "?"), so
they read as thermodynamically undefined. This script closes that gap: for every
reaction whose canonical deltag is missing, it promotes the best available
per-source estimate into deltag / deltagerr and adopts that estimate's own
direction operator as the reversibility.

It is purely a re-aggregation of values already in the database — no new
estimation, no external dependencies. Reactions that already have a canonical
deltag are left untouched.

Selection policy: prefer the mechanistic, experimentally-anchored tier
(eQuilibrator, then Group Contribution) over the machine-learning tier
(dGPredictor-ModelSEED, dGPredictor); WITHIN the chosen tier take the estimate
with the smallest reported uncertainty (ties broken by the listed order). The
within-tier lowest-error rule matters: it stops a wildly-uncertain ML outlier
(e.g. -100 +/- 71 kcal/mol) from being promoted over a tight prediction
(-8.6 +/- 0.04) from the same tier. The lower-confidence de-novo eQuilibrator
source is intentionally excluded. Edit TIERS to change the policy.
"""
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions
from Estimate_Reaction_Reversibility import reversibility_from_energy

SENTINEL = 10000000
MAX_ABS_DG = 1000.0   # kcal/mol; reject implausible magnitudes
MAX_ERR = 100.0       # kcal/mol; reject uselessly-uncertain estimates

# Quality tiers, best first; within a tier the lowest-error usable estimate wins
# (the per-tier list order is only a tie-breaker).
TIERS = [
    ["eQuilibrator", "Group contribution"],          # mechanistic / measurement-anchored
    ["dGPredictor-ModelSEED", "dGPredictor"],         # machine-learning predictors
]


def is_missing(dg):
    if dg is None:
        return True
    try:
        return abs(float(dg)) >= SENTINEL
    except (TypeError, ValueError):
        return True


def usable(entry):
    """A per-source [energy, error, operator] is usable for promotion iff the
    energy is a finite, non-sentinel, plausible number with a bounded error."""
    if not isinstance(entry, list) or len(entry) < 2:
        return False
    dg, err = entry[0], entry[1]
    if not isinstance(dg, (int, float)) or isinstance(dg, bool):
        return False
    if dg != dg or abs(dg) == float("inf") or abs(dg) >= SENTINEL or abs(dg) > MAX_ABS_DG:
        return False
    if not isinstance(err, (int, float)) or isinstance(err, bool) or err != err:
        err = 0.0
    if abs(err) >= SENTINEL or abs(err) > MAX_ERR:
        return False
    return True


def main():
    helper = Reactions()
    reactions = helper.loadReactions()

    promoted = 0
    by_source = {}
    skipped_no_source = 0
    examined = 0

    for rxn in sorted(reactions.keys()):
        robj = reactions[rxn]
        if robj.get('status') == 'EMPTY':
            continue
        if not is_missing(robj.get('deltag')):
            continue  # canonical value already present — never overwrite it
        examined += 1

        thermo = robj.get('thermodynamics')
        if not isinstance(thermo, dict):
            skipped_no_source += 1
            continue

        chosen = None
        for tier in TIERS:
            in_tier = [s for s in tier if s in thermo and usable(thermo[s])]
            if in_tier:
                # smallest reported uncertainty wins; tie-break by tier order
                chosen = min(in_tier, key=lambda s: (abs(thermo[s][1]), tier.index(s)))
                break
        if chosen is None:
            skipped_no_source += 1
            continue

        entry = thermo[chosen]
        dg = float("{0:.2f}".format(float(entry[0])))
        err = float("{0:.2f}".format(float(entry[1])))

        # Adopt the source's own direction operator (computed by the same
        # heuristic as the canonical reversibility); recompute as a fallback.
        op = entry[2] if (len(entry) >= 3 and entry[2] in (">", "<", "=")) else None
        if op is None:
            op = reversibility_from_energy(robj, dg, err)

        robj['deltag'] = dg
        robj['deltagerr'] = err
        robj['reversibility'] = op
        promoted += 1
        by_source[chosen] = by_source.get(chosen, 0) + 1

    print("Reactions with a missing canonical deltag examined: " + str(examined))
    print("Promoted to a canonical deltag: " + str(promoted))
    for src in [s for tier in TIERS for s in tier]:
        if src in by_source:
            print("  from " + src + ": " + str(by_source[src]))
    print("Left undefined (no usable per-source estimate): " + str(skipped_no_source))
    print("Saving reactions")
    helper.saveReactions(reactions)


if __name__ == "__main__":
    main()
