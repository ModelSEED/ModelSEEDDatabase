#!/usr/bin/env python3
"""Audit atom mappings against the one class where chemistry gives a hard answer.

Run from the repository root with an interpreter that has rdkit, e.g.
    /scratch/seaver/micromamba/envs/equilibrator/bin/python \
        Scripts/Structures/audit_decarboxylation_mappings.py

The test
========

A carbon that becomes CO2 must have carried two oxygens. A carbon with none
cannot possibly become CO2. For every reaction producing CO2, find the
substrate carbon mapped to CO2's carbon and count that carbon's oxygen
neighbours. Zero oxygens is a proven-wrong mapping, not a heuristic.

Single-oxygen carbons are NOT counted as errors -- oxidative decarboxylation
can legitimately recruit an oxygen from water -- so the impossible count is a
floor on the true error rate, not an estimate of it.

What it currently reports (2026-09-14, after the chain-salvage withdrawal)
=========================================================================

  745 reactions produce CO2 with a carbon mapped to it
  488 (57.1%) correct, 181 (21.2%) possible, 185 (21.7%) IMPOSSIBLE
  185 impossible rows across 174 distinct reactions

Two conclusions, both load-bearing:

1. InChI order is the right reading of the `E#N` indices (it fits roughly twice
   as well as SMILES order), confirming AtomMappings/README.md.

2. The errors belong to RDT, not to this repository. Verified by rerunning RDT
   on a 13-reaction panel and decoding its native ECBLAST atom-atom-mapping
   numbers directly: RDT's own output and the shipped rows agreed on 376 of 376
   atom pairs. rxn00168 (pyruvate decarboxylase) is the cleanest illustration --
   RDT emits
       [O:5]=[C:4]([O-:6])[C:2](=[O:3])[CH3:1]>>[C:1](=[O:3])=[O:6].[CH3:2][CH:4]=[O:5]
   i.e. it maps pyruvate's METHYL carbon to CO2 instead of the carboxylate.
   Deterministic, flag-independent, and green under upstream's own EC 4.1.1.1
   test, which only asserts that the bond-change fingerprint is non-empty.

   This audit previously reported 187 rows / 176 reactions. Two of those were
   artifacts of this repository's own chain-splitting defect and disappeared
   when that was withdrawn; the remaining 174 are RDT's.

Caveat on the confidence tag
============================

Within this decarboxylation set `clean` fails at about the same rate as
`salvaged`, so the tag does not predict correctness HERE. Do not generalise
that: across the whole set the tag is strongly predictive (roughly 29% of
`clean` reactions flagged by structural validators versus 92% of `salvaged`).
This class is enriched for RDT's known weak spot.
"""
import csv, json, glob, collections
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

STRUCT = 'Biochemistry/Structures/Unique_ModelSEED_Structures.txt'
MAP    = 'Biochemistry/Structures/AtomMappings/all_mapping_no_problem.txt'
CO2    = 'cpd00011'

smiles, inchi = {}, {}
with open(STRUCT) as fh:
    for row in csv.reader(fh, delimiter='\t'):
        if len(row) < 6: continue
        if row[1] == 'SMILE': smiles[row[0]] = row[5]
        elif row[1] == 'InChI': inchi[row[0]] = row[5]

def carbon_oxygen_counts(cid, order):
    """[n_oxygen_neighbours] indexed by per-element carbon position."""
    if order == 'smiles':
        m = Chem.MolFromSmiles(smiles.get(cid, '')) if cid in smiles else None
    else:
        m = Chem.MolFromInchi(inchi.get(cid, '')) if cid in inchi else None
    if m is None: return None
    out = []
    for a in m.GetAtoms():
        if a.GetSymbol() != 'C': continue
        out.append(sum(1 for nb in a.GetNeighbors() if nb.GetSymbol() == 'O'))
    return out

# which substrate carbon maps to CO2's carbon, per reaction
tgt = collections.defaultdict(list)
for line in open(MAP):
    rx, rest = line.split(' ', 1)
    a, b = rest.strip().split('=')
    sa, ea = a.split(':'); sb, eb = b.split(':')
    if sb != CO2 or not eb.startswith('C#') or not ea.startswith('C#'): continue
    if sa == CO2: continue
    tgt[rx].append((sa, int(ea[2:])))

res = {o: collections.Counter() for o in ('inchi', 'smiles')}
bad = {o: [] for o in ('inchi', 'smiles')}
for rx, lst in tgt.items():
    for cid, idx in lst:
        for order in ('inchi', 'smiles'):
            counts = carbon_oxygen_counts(cid, order)
            if counts is None or idx > len(counts):
                res[order]['unparsable'] += 1; continue
            n_ox = counts[idx-1]
            if n_ox >= 2:   res[order]['carboxyl (2+ O) - correct'] += 1
            elif n_ox == 1: res[order]['1 oxygen - possible'] += 1
            else:
                res[order]['0 oxygen - IMPOSSIBLE'] += 1
                bad[order].append((rx, cid, idx))

print(f"  reactions producing CO2 with a mapped carbon: {len(tgt):,}\n")
for order in ('inchi', 'smiles'):
    print(f"  --- indices read as {order.upper()} order ---")
    t = sum(res[order].values())
    for k, v in sorted(res[order].items(), key=lambda z: -z[1]):
        print(f"      {k:<28} {v:>6,}  ({100*v/t:.1f}%)")
    print(f"      examples of impossible: {[b[0] for b in bad[order][:6]]}\n")

conf = {r.split('\t')[0]: r.split('\t')[1].strip()
        for r in open('Biochemistry/Structures/AtomMappings/rxns_confidence.tsv')}
rxbad = {b[0] for b in bad['inchi']}
allrx = set(tgt)
print("  === does atom_mapping_confidence predict the error? ===")
for tag in ('clean', 'salvaged'):
    tot = sum(1 for r in allrx if conf.get(r) == tag)
    err = sum(1 for r in rxbad if conf.get(r) == tag)
    if tot: print(f"    {tag:<9} {err:>4,} impossible of {tot:>4,}  ({100*err/tot:.1f}%)")
print(f"\n  rxn00168 confidence: {conf.get('rxn00168')}   in the impossible set: {'rxn00168' in rxbad}")
