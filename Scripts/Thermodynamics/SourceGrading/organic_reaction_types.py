#!/usr/bin/env python3
"""Classify a ModelSEED reaction by its ORGANIC-CHEMISTRY transformation type.

EC class answers "what enzyme family is this" and is only a loose proxy for
"what bond changes". EC 1 lumps hydride transfer, O2 insertion, disulfide
exchange and quinone reduction together, and those are four different chemical
problems for a group-contribution estimator. This module assigns each reaction
one transformation class based on what actually changes: which cofactor couple
carries the electrons, whether O2 is incorporated or reduced to water, which
bond is made or broken, and which group is moved.

Classification is a PRIORITY CASCADE -- first match wins -- ordered
most-chemically-specific first, so a NAD-linked quinone reductase is called
quinoid redox rather than hydride transfer. Each reaction gets exactly one
label, which is what makes the class shares interpretable as a partition.

Inputs per reaction: the ModelSEED stoichiometry, the compound records, the EC
numbers, and the RDKit functional-group deltas already computed into
reaction_features.tsv (``d_<group>`` = net count of that group created).

``classify(row, stoich, cpds)`` -> (class, subclass).
"""
from __future__ import annotations

import re

# --- cofactor / carrier identities -----------------------------------------
NADS = {"cpd00003", "cpd00004", "cpd00005", "cpd00006"}          # NAD(P)(H)
FLAVINS = {"cpd00015", "cpd00982", "cpd00050", "cpd01270"}        # FAD(H2), FMN(H2)
FERREDOXIN = {"cpd11620", "cpd11621"}
CYTOCHROME = {"cpd00109", "cpd00110"}
O2, H2O2, H2O, CO2 = "cpd00007", "cpd00025", "cpd00001", "cpd00011"
ATP, ADP, AMP, PI, PPI = "cpd00002", "cpd00008", "cpd00018", "cpd00009", "cpd00012"
SAM, SAH = "cpd00017", "cpd00019"
COA = "cpd00010"
GLU, GLN, ASP, ASN = "cpd00023", "cpd00053", "cpd00041", "cpd00132"
THF_POOL = {"cpd00087", "cpd00345", "cpd00443", "cpd00125"}

# Quinone/quinol pairs are matched on name because ModelSEED carries dozens of
# prenyl-chain-length variants (ubiquinone-6..10, menaquinone-6..9, plasto-,
# phyllo-) that share no id prefix.
QUINONE_RE = re.compile(
    r"(ubiquino|menaquino|plastoquino|phylloquino|naphthoquino|naphtoquin|"
    r"ubidecarenone|ubiquinol|menaquinol|plastoquinol|demethylmenaquino|"
    r"vitamin k)", re.I)
# Catechol-type dihydroxyaromatics: the same 2e-/2H+ quinoid couple on a ring
CATECHOLIC_RE = re.compile(
    r"(catechol|protocatechu|homogentisate|dihydroxybenzo|dihydroxyphenyl|"
    r"hydroquinone|dopamine|caffeate|caffeoyl|quinol)", re.I)


def _names(stoich: dict, cpds: dict) -> list[str]:
    return [str(cpds.get(c, {}).get("name", "")) for c in stoich]


def classify(row, stoich: dict, cpds: dict) -> tuple[str, str]:
    """Return (organic transformation class, subclass) for one reaction."""
    ids = set(stoich)
    ec = str(row.get("ec", "") or "")
    ecs = [e for e in ec.split(";") if e]
    names = _names(stoich, cpds)
    n_o2 = -float(stoich.get(O2, 0.0))          # >0 means O2 consumed

    def d(g):                                    # net change in a group count
        try:
            return float(row.get(f"d_{g}", 0.0) or 0.0)
        except (TypeError, ValueError):
            return 0.0

    def ec_is(*prefixes):
        return any(e.startswith(p) for e in ecs for p in prefixes)

    has_quinone = any(QUINONE_RE.search(n) for n in names)
    has_catecholic = any(CATECHOLIC_RE.search(n) for n in names)

    # 1. Quinoid / aromatic two-electron redox. Kept first because it is the
    #    class the two estimators disagree on most and it hides inside EC 1.
    if has_quinone:
        return ("Redox: quinone / quinol", "prenyl-quinone carrier")
    if has_catecholic and (d("hydroxyl") != 0 or d("ketone") != 0 or n_o2 > 0):
        return ("Redox: quinone / quinol", "catecholic (dihydroxyarene)")

    # 2. O2 chemistry, split by fate of the oxygen atoms.
    if n_o2 > 0:
        if ec_is("1.13"):
            return ("Oxygenation: O2 incorporated", "dioxygenase (both O atoms)")
        if ec_is("1.14"):
            return ("Oxygenation: O2 incorporated", "monooxygenase / hydroxylase")
        if H2O2 in ids or ec_is("1.1.3", "1.4.3", "1.5.3", "1.2.3", "1.7.3", "1.8.3"):
            return ("Redox: O2 as terminal acceptor", "oxidase -> H2O2 / H2O")
        return ("Oxygenation: O2 incorporated", "unclassified oxygen transfer")

    # 3. Heteroatom redox -- sulfur and nitrogen behave unlike carbon.
    if d("disulfide") != 0 or d("thiol") != 0:
        return ("Redox: heteroatom (S)", "thiol / disulfide")
    if d("nitro") != 0 or d("nitrile") != 0:
        return ("Redox: heteroatom (N)", "nitro / nitrile")

    # 4. Carbon redox by the carrier that takes the electrons.
    if ids & NADS:
        return ("Redox: carbon, hydride transfer", "NAD(P)-linked")
    if ids & FLAVINS:
        return ("Redox: carbon, hydride transfer", "flavin-linked")
    if ids & (FERREDOXIN | CYTOCHROME):
        return ("Redox: carbon, hydride transfer", "ferredoxin / cytochrome")

    # 5. Group transfer -- a conserved moiety moves between two carriers.
    if d("phosphoanhydride") != 0 or (ids & {ATP, ADP, AMP} and (PI in ids or PPI in ids)):
        return ("Group transfer: phosphoryl", "kinase / phosphatase / anhydride")
    if ids & {SAM, SAH} or (ids & THF_POOL and ec_is("2.1.1")):
        return ("Group transfer: methyl / C1", "SAM or folate C1")
    if d("thioester") != 0 or COA in ids:
        return ("Group transfer: acyl (thioester)", "CoA / ACP acyl")
    if d("acetal") != 0 and ec_is("2.4", "3.2"):
        return ("Group transfer: glycosyl", "glycosyl / glycosidase")
    if ids & {GLU, GLN, ASP, ASN} and (d("primary_amine") != 0 or d("imine") != 0
                                       or ec_is("2.6.1", "6.3.5")):
        return ("Group transfer: amino / amido", "transaminase / amidotransferase")

    # 6. C-C bond changes.
    if CO2 in ids or ec_is("4.1.1", "6.4.1"):
        return ("C-C bond: carboxylation / decarboxylation", "CO2 added or lost")
    if ec_is("4.1.2", "2.3.3", "4.1.3"):
        return ("C-C bond: aldol / Claisen", "condensation or retro-cleavage")

    # 7. Additions and eliminations across a C=C / C=O.
    if ec_is("4.2.1"):
        return ("Addition / elimination", "hydratase / dehydratase")
    if ec_is("4.2", "4.3", "4.4", "4.5", "4.6", "4.99"):
        return ("Addition / elimination", "other lyase")

    # 8. Hydrolysis -- water cleaves a bond, no redox, no carrier.
    if ec_is("3.") or (H2O in ids and float(stoich.get(H2O, 0.0)) < 0
                       and (d("ester") != 0 or d("amide") != 0 or d("acetal") != 0)):
        sub = ("ester" if d("ester") != 0 else
               "amide / peptide" if d("amide") != 0 else
               "glycosidic" if d("acetal") != 0 else "other hydrolase")
        return ("Hydrolysis", sub)

    # 9. Skeletal rearrangement with no net redox or group exchange.
    if ec_is("5."):
        sub = ("racemase / epimerase" if ec_is("5.1") else
               "cis-trans isomerase" if ec_is("5.2") else
               "intramolecular oxidoreductase" if ec_is("5.3") else
               "mutase" if ec_is("5.4") else "other isomerase")
        return ("Isomerisation / rearrangement", sub)

    # 10. ATP-driven bond formation not already caught as phosphoryl transfer.
    if ec_is("6."):
        return ("Ligation (ATP-driven bond formation)", "ligase")

    return ("Other / unclassified", "-")


CLASS_ORDER = [
    "Redox: quinone / quinol",
    "Oxygenation: O2 incorporated",
    "Redox: O2 as terminal acceptor",
    "Redox: heteroatom (S)",
    "Redox: heteroatom (N)",
    "Redox: carbon, hydride transfer",
    "Group transfer: phosphoryl",
    "Group transfer: methyl / C1",
    "Group transfer: acyl (thioester)",
    "Group transfer: glycosyl",
    "Group transfer: amino / amido",
    "C-C bond: carboxylation / decarboxylation",
    "C-C bond: aldol / Claisen",
    "Addition / elimination",
    "Hydrolysis",
    "Isomerisation / rearrangement",
    "Ligation (ATP-driven bond formation)",
    "Other / unclassified",
]
