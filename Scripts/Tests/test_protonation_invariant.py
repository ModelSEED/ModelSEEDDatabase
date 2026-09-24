#!/usr/bin/env python3
"""Regression test: no protonation row ships that a protonation cannot produce.

Guards three things:

1. **The invariant itself** -- ``compare`` must call a proton moved
   (dH == dcharge) fine, and call the two ways this repository has broken it
   by name: hydrogens gained with no charge change (InChI disconnected a metal
   cluster and Marvin protonated the freed ligands: Fe4S4 -> H8Fe4S4), and
   charge changed with no hydrogen change.

2. **InChI layer derivation** -- ``inchi_layers`` must return what the string
   declares, on the cases where a molecule parser did not: acetate (/p-1),
   chlorate (/q-1 on a hypervalent chlorine that RDKit rejects), a shattered
   iron-sulfur cluster, triphenyltin chloride (five InChI components), and a
   magnesium porphyrin with both /q and /p layers.

3. **The shipped bundles** -- every consumed protonation bundle must pass the
   ``row`` check with zero strict failures. This is the check that would have
   caught 37 impossible formulas before they reached the compound records.

Exit code 0 on all-pass, 1 on any failure. Suitable for CI.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "Structures"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "Libs", "Python"))

import Validate_Protonations as V                          # noqa: E402
from Print_Structure_Formula_Charge import inchi_layers    # noqa: E402

# Bundles that feed the compound records. Both are globbed by
# BiochemPy.loadStructures, so both must be clean.
BUNDLES = ("marvin_26.1_ph7", "marvin_23.4_ph7")

failures = []


def check(name, cond, detail=""):
    if not cond:
        failures.append(f"{name}: {detail}")


# 1. the invariant
check("proton removed", V.compare("C2H4O2", "0", "C2H3O2", "-1")[0] == "ok")
check("two protons removed", V.compare("H3O4P", "0", "HO4P", "-2")[0] == "ok")
check("unchanged", V.compare("Fe4S4", "2", "Fe4S4", "2")[0] == "ok")
check("Fe4S4 shattered", V.compare("Fe4S4", "0", "H8Fe4S4", "0")[0] == "H gained, no charge change",
      V.compare("Fe4S4", "0", "H8Fe4S4", "0"))
check("triphenyltin", V.compare("C18H15ClSn", "0", "C18H18ClSn", "0")[0] == "H gained, no charge change")
check("elemental S", V.compare("S", "0", "H2S", "0")[0] == "H gained, no charge change")
check("charge without H", V.compare("ClO2", "0", "ClO2", "-1")[0] == "charge changed, no H change")
check("H lost without charge", V.compare("HClO3", "0", "ClO3", "0")[0] == "H lost, no charge change")
check("heavy atoms", V.compare("Fe4S4", "0", "Fe4S3", "0")[0].startswith("heavy atoms changed"))
check("wildcard is INFO", V.compare("H4NR", "1", "H3NR", "1")[0] == "wildcard")
check("missing source is INFO", V.compare("", "", "H2S", "0")[0] == "unchecked")
check("INFO kinds do not fail", all(k in V.INFO_KINDS for k in ("wildcard", "unchecked")))

# 2. InChI layers
cases = {
    "acetate":       ("InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)/p-1", ("C2H3O2", -1)),
    "chlorate":      ("InChI=1S/ClO3/c2-1(3)4/q-1", ("ClO3", -1)),
    "chlorite":      ("InChI=1S/ClO2/c2-1-3/q-1", ("ClO2", -1)),
    "Fe4S4 +1":      ("InChI=1S/4Fe.4S/q;;;+1;;;;", ("Fe4S4", 1)),
    "Fe4S4 shattered": ("InChI=1S/4Fe.4H2S/h;;;;4*1H2", ("H8Fe4S4", 0)),
    "triphenyltin":  ("InChI=1S/3C6H5.ClH.Sn/c3*1-2-4-6-5-3-1;;/h3*1-5H;1H;/q;;;;+1/p-1", ("C18H15ClSn", 0)),
    "Mg porphyrin":  ("InChI=1S/C35H35N4O4.Mg/c1-8-22-18(3)26-14-27-20(5)24(10-12-34(40)41)32(38-27)"
                      "17-33-25(11-13-35(42)43-7)21(6)29(39-33)16-31-23(9-2)19(4)28(37-31)15-30(22)36-26;"
                      "/h8-9,14-17H,1-2,10-13H2,3-7H3,(H2-,36,37,38,39,40,41);/q-1;+2/p-1",
                      ("C35H34MgN4O4", 0)),
    "multiplied q":  ("InChI=1S/2ClO3.Mg/c2*2-1(3)4;/q2*-1;+2", ("Cl2MgO6", 0)),
}
for name, (inchi, want) in cases.items():
    got = inchi_layers(inchi)
    check(f"layers {name}", got == want, f"got {got}, want {want}")
check("non-standard InChI falls back", inchi_layers("InChI=1/C2H4O2/c1-2(3)4/h1H3,(H,3,4)") == (None, None))
check("fixed-H layer falls back", inchi_layers("InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)/f/h3H") == (None, None))

# 3. the shipped bundles
allowed = V.load_allowlist()
for bundle in BUNDLES:
    summary, findings = V.run(bundle, list(V.SOURCES), want_reactions=False)
    strict = [f for f in findings if f["check"] == "row" and f["kind"] not in V.INFO_KINDS
              and (f["source"], f["external_id"], "row") not in allowed]
    layers = [f for f in findings if f["check"] == "layers"]
    keys = [f for f in findings if f["check"] == "key"]
    check(f"{bundle}: every InChIKey row is the key of its InChI row", not keys,
          f"{len(keys)} mismatches, e.g. " + "; ".join(
              f"{f['source']}:{f['external_id']}" for f in keys[:3]))
    check(f"{bundle}: rows examined", sum(v for (s, c, t), v in summary.items() if c == "row") > 90000,
          f"only {sum(v for (s, c, t), v in summary.items() if c == 'row')} rows examined")
    check(f"{bundle}: no impossible protonation rows", not strict,
          f"{len(strict)} strict row failures, e.g. " + "; ".join(
              f"{f['source']}:{f['external_id']} {f['type']} {f['from_formula']}/{f['from_charge']}"
              f"->{f['to_formula']}/{f['to_charge']}" for f in strict[:3]))
    check("inchi.tsv agrees with its own InChI layers", not layers,
          f"{len(layers)} rows disagree, e.g. " + "; ".join(
              f"{f['source']}:{f['external_id']}" for f in layers[:3]))

if failures:
    print("FAIL")
    for f in failures:
        print("  -", f)
    sys.exit(1)
print(f"OK: invariant, InChI layers, and {len(BUNDLES)} bundles clean")
