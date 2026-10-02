#!/usr/bin/env python
"""Numbers for reviewer 1 comment 2 (transport) and reviewer 2 comment 2 (LLMs).

Both questions are about what a claim in the paper is actually worth, and both
are answerable from the shipped Biochemistry/reaction_*.json alone -- no
eQuilibrator working tree needed, unlike the other scripts in this directory.

Transport: the paper says transport reactions are "scored from stoichiometry
and energy alone", with no membrane potential and no pH gradient. The reviewer
asks what that costs. Quantify it: how many transport reactions each source
commits on, and what grade they end up carrying, against the non-transport
remainder.

LLMs: the reviewer cannot see what the ensemble is for. The defensible answer
is coverage -- the reactions the thermodynamic rules leave undirected. Measure
that, plus agreement with eQuilibrator where both commit AND against the
measured energies in opentecr_comparison.csv, so the section can state a role
instead of a motivation.

POPULATION. Every count here includes obsolete records, because that is what the
manuscript counts: M12's "30,157 ... and 25,855 receive none" sums to 56,012,
and the grade totals reproduce Supplementary Table S1 exactly on this basis
(33,099 graded; 3,434/18,388/11,277; 806 anchors). Filtering obsolete rows out
gives 27,355/2,486/15,520/9,349/365 and silently disagrees with the paper --
which is a filter artefact, not a rebuild regression. Pass --live to see that
basis; do not quote it in the manuscript without converting every other number
with it.

Run:  ~/Documents/py_venv/bin/python review_transport_and_llm.py
"""
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PREDICTORS = ("eQuilibrator", "dGPredictor", "Group contribution")
COMMITTED = (">", "<")          # an actual direction call
STATED = (">", "<", "=")        # a direction or an explicit reversibility call


LIVE_ONLY = "--live" in sys.argv


def load():
    for shard in sorted(ROOT.glob("Biochemistry/reaction_[0-9][0-9].json")):
        for rxn in json.load(open(shard)):
            if LIVE_ONLY and rxn.get("is_obsolete") == 1:
                continue
            yield rxn


def direction(rxn, source):
    """Third slot of the per-source thermodynamics triple, or None."""
    t = (rxn.get("thermodynamics") or {}).get(source)
    if not t or len(t) < 3:
        return None
    return t[2] or None


def has_energy(rxn, source):
    t = (rxn.get("thermodynamics") or {}).get(source)
    return bool(t) and t[0] not in ("", None)


def pct(n, d):
    return f"{100 * n / d:.1f}%" if d else "  n/a"


def main():
    rxns = list(load())
    print(f"population: {'live only (--live)' if LIVE_ONLY else 'ALL records, the basis the manuscript uses'}")
    tr = [r for r in rxns if r.get("is_transport") == 1]
    nt = [r for r in rxns if r.get("is_transport") != 1]

    print("=" * 78)
    print("TRANSPORT REACTIONS  (reviewer 1, comment 2)")
    print("=" * 78)
    print(f"reactions               {len(rxns):>8,}")
    print(f"  transport             {len(tr):>8,}  ({pct(len(tr), len(rxns))})")
    print(f"  non-transport         {len(nt):>8,}")

    print(f"\n{'':<22}{'transport':>22}{'non-transport':>22}")
    print(f"{'':<22}{'n':>8}{'share':>8}{'':>6}{'n':>8}{'share':>8}")
    print("-" * 66)

    def row(label, ftr, fnt):
        a = sum(1 for r in tr if ftr(r))
        b = sum(1 for r in nt if fnt(r))
        print(f"{label:<22}{a:>8,}{pct(a, len(tr)):>8}{'':>6}{b:>8,}{pct(b, len(nt)):>8}")

    for s in PREDICTORS:
        row(f"{s[:20]} energy", lambda r, s=s: has_energy(r, s),
            lambda r, s=s: has_energy(r, s))
    print()
    for s in PREDICTORS:
        row(f"{s[:20]} dir.", lambda r, s=s: direction(r, s) in COMMITTED,
            lambda r, s=s: direction(r, s) in COMMITTED)
    print()
    row("any source states", lambda r: any(direction(r, s) in STATED for s in PREDICTORS),
        lambda r: any(direction(r, s) in STATED for s in PREDICTORS))
    row("graded at all", lambda r: bool(r.get("thermo-evidence")),
        lambda r: bool(r.get("thermo-evidence")))

    print("\ngrade split")
    for label, group in (("transport", tr), ("non-transport", nt)):
        g = Counter((r.get("thermo-evidence") or {}).get("grade") for r in group)
        tot = sum(v for k, v in g.items() if k)
        line = "  ".join(f"{k} {g.get(k, 0):,} ({pct(g.get(k, 0), tot)})"
                         for k in ("gold", "silver", "bronze"))
        print(f"  {label:<16}{line}   graded n={tot:,}")

    print("\nwhat the transport grades actually rest on")
    both = sum(1 for r in tr if _spans(r))
    print(f"  stoichiometry lists two compartments  {both:>7,}  ({pct(both, len(tr))})")
    uni = [r for r in tr if _uniport(r)]
    ug = Counter((r.get("thermo-evidence") or {}).get("grade") for r in uni)
    print(f"  uniport-like, no net chemistry        {len(uni):>7,}  ({pct(len(uni), len(tr))})")
    print(f"    their grades: " + ", ".join(
        f"{k or 'ungraded'} {v:,}" for k, v in ug.most_common()))

    # The gold transport reactions are the ones worth explaining: an energy-only
    # score cannot see a membrane, so a gold grade here is confidence in the
    # COUPLED CHEMISTRY, not in the translocation.
    gold = [r for r in tr if (r.get("thermo-evidence") or {}).get("grade") == "gold"]
    atp = sum(1 for r in gold if {"cpd00002", "cpd00008"} <= _cpds(r))
    hplus = sum(1 for r in gold if _crosses(r, "cpd00067"))
    print(f"\n  gold-graded transport                 {len(gold):>7,}  "
          f"({pct(len(gold), len(tr))} of transport; "
          f"{pct(sum(1 for r in nt if (r.get('thermo-evidence') or {}).get('grade') == 'gold'), len(nt))} of non-transport)")
    print(f"    carry ATP + ADP (coupled hydrolysis) {atp:>7,}  ({pct(atp, len(gold))})")
    print(f"    translocate a proton themselves      {hplus:>7,}  ({pct(hplus, len(gold))})")
    print(f"    deciding source: " + ", ".join(
        f"{k} {v:,}" for k, v in Counter(
            (r.get("thermo-evidence") or {}).get("source") for r in gold).most_common()))

    print("\n" + "=" * 78)
    print("LLM ENSEMBLE  (reviewer 2, comment 2)")
    print("=" * 78)
    llm = {r["id"]: direction(r, "LLMs") for r in rxns
           if direction(r, "LLMs") is not None}
    print(f"reactions carrying an LLM call      {len(llm):>8,}  "
          f"({pct(len(llm), len(rxns))})")

    thermo_dir = {r["id"] for r in rxns
                  if any(direction(r, s) in STATED for s in PREDICTORS)}
    none_dir = {r["id"] for r in rxns} - thermo_dir
    gain = none_dir & set(llm)
    print(f"reactions NO predictor directs      {len(none_dir):>8,}")
    print(f"  ... of which the LLMs do direct   {len(gain):>8,}  "
          f"({pct(len(gain), len(none_dir))})  <- the coverage argument")

    # agreement with eQuilibrator where both commit to an irreversible call
    eq = {r["id"]: direction(r, "eQuilibrator") for r in rxns
          if direction(r, "eQuilibrator") in COMMITTED}
    both_commit = [i for i in eq if llm.get(i) in COMMITTED]
    agree = sum(1 for i in both_commit if eq[i] == llm[i])
    print(f"\nboth eQuilibrator and LLMs commit   {len(both_commit):>8,}")
    print(f"  agree on the direction            {agree:>8,}  ({pct(agree, len(both_commit))})")

    # agreement restricted to the measured anchors (gold, source 'measured')
    anch = [r for r in rxns
            if (r.get("thermo-evidence") or {}).get("assessment") == "measured"]
    print(f"\nmeasured anchors                    {len(anch):>8,}")
    ac = [r for r in anch if llm.get(r["id"]) in COMMITTED]
    print(f"  LLMs commit on                    {len(ac):>8,}")
    # Denominator must be anchors where BOTH commit -- scoring the LLM against
    # an eQuilibrator call that was never made reads as disagreement when it is
    # simply eQuilibrator abstaining.
    bc = [r for r in ac if direction(r, "eQuilibrator") in COMMITTED]
    ok = sum(1 for r in bc if llm[r["id"]] == direction(r, "eQuilibrator"))
    print(f"  ... and eQuilibrator also commits {len(bc):>8,}")
    print(f"      the two agree                 {ok:>8,}  ({pct(ok, len(bc))})")

    print("\nLLM call distribution")
    for k, v in Counter(llm.values()).most_common():
        print(f"  {k!r:<6}{v:>8,}  ({pct(v, len(llm))})")

    _vs_measured(llm)

    print("\nLLM coverage on transport specifically")
    trl = sum(1 for r in tr if r["id"] in llm)
    print(f"  transport with an LLM call        {trl:>8,}  ({pct(trl, len(tr))})")


def _vs_measured(llm):
    """Score the ensemble against MEASUREMENT, not against another prediction.

    opentecr_comparison.csv ships opentecr_dG_kJ already written in the
    ModelSEED equation's orientation -- confirmed by the shipped eQuilibrator
    energy matching its sign on 183 of 185 anchors. Do NOT apply
    ms_orientation_vs_canonical as a flip: that drops eQuilibrator to 45% and
    is the wrong reading of the column.
    """
    path = ROOT / "Biochemistry/Thermodynamics/SourceGrading/opentecr_comparison.csv"
    if not path.exists():
        return
    rows = list(csv.DictReader(open(path)))
    by_id = {}
    for shard in sorted(ROOT.glob("Biochemistry/reaction_[0-9][0-9].json")):
        for r in json.load(open(shard)):
            by_id[r["id"]] = r
    print("\ndirection vs MEASURED energy (openTECR), by decision margin")
    print(f"  {'margin':<12}{'source':<16}{'n':>6}{'agree':>7}{'rate':>8}{'chance':>8}{'kappa':>7}")
    for margin in (0.0, 5.0, 11.7):
        for src in ("eQuilibrator", "LLMs"):
            pairs = []
            for row in rows:
                rxn = by_id.get(row["modelseed_rxn"])
                if not rxn:
                    continue
                try:
                    g = float(row["opentecr_dG_kJ"])
                except (TypeError, ValueError):
                    continue
                m = ">" if g < -margin else ("<" if g > margin else "=")
                c = direction(rxn, src)
                if m in COMMITTED and c in COMMITTED:
                    pairs.append((m, c))
            if not pairs:
                continue
            n = len(pairs)
            agree = sum(1 for a, b in pairs if a == b)
            pm, pc = Counter(a for a, _ in pairs), Counter(b for _, b in pairs)
            chance = sum((pm[d] / n) * (pc[d] / n) for d in COMMITTED)
            kappa = (agree / n - chance) / (1 - chance) if chance < 1 else float("nan")
            print(f"  {margin:<12}{src:<16}{n:>6}{agree:>7}{100*agree/n:>7.1f}%"
                  f"{100*chance:>7.1f}%{kappa:>7.2f}")
    # the confusion matrix is the finding: the errors are all on one side
    cm = Counter()
    for row in rows:
        rxn = by_id.get(row["modelseed_rxn"])
        if not rxn:
            continue
        try:
            g = float(row["opentecr_dG_kJ"])
        except (TypeError, ValueError):
            continue
        m = ">" if g < -11.7 else ("<" if g > 11.7 else "=")
        c = direction(rxn, "LLMs")
        if m in COMMITTED and c in COMMITTED:
            cm[(m, c)] += 1
    print("  confusion at 11.7 kJ (measured -> LLM): "
          + ", ".join(f"{a}->{b} {n}" for (a, b), n in sorted(cm.items())))


def _spans(rxn):
    """True when the stoichiometry names more than one compartment index."""
    return len({p["compartment"] for p in rxn.get("stoichiometry") or []}) > 1


def _uniport(rxn):
    """True when every species in the reaction merely changes compartment -- the
    case whose net chemical stoichiometry is empty, so an energy-only score has
    no signal at all. Checked, and these are NOT the gold ones: the scheme
    already sends them to bronze or leaves them ungraded."""
    by_cpd = defaultdict(set)
    for p in rxn.get("stoichiometry") or []:
        by_cpd[p["compound"]].add(p["compartment"])
    moved = {c for c, cp in by_cpd.items() if len(cp) > 1}
    return bool(moved) and moved == set(by_cpd)


def _cpds(rxn):
    return {p["compound"] for p in rxn.get("stoichiometry") or []}


def _crosses(rxn, cpd):
    """True when this compound appears in more than one compartment."""
    return len({p["compartment"] for p in rxn.get("stoichiometry") or []
                if p["compound"] == cpd}) > 1


if __name__ == "__main__":
    sys.exit(main())
