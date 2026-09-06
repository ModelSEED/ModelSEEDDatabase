# AtomMappings — open considerations

Working-notes companion to `README.md` and
`Papers/NAR_Update_2026/drafts/methods_atom_mapping.md`. Tracks the
issues and design questions that were surfaced but deliberately not
resolved during the 2026 atom-mapping arc. Live document — update as
new items surface or old ones close out.

## Things that work today

- RDT-produced mapping rows are captured for **32,877 / 56,012** reactions
  (59% full DB; 80.9% priority scope; 93.1% Athaliana_TAIR10 model),
  ingested onto reaction JSONs as a single `atom_mapping` dict with
  `data`, `confidence`, and `has_symmetry_groups` fields.
- Row-level filter (`Rebuild_AtomMappings_from_raw.py`) recovers atom
  pairs from RDT's run-on chains / dangling orphans / cross-element
  rows that Sebastian's original whole-reaction reject-on-any-problem
  filter dropped wholesale (contributed upstream).
- Confidence tagging distinguishes `clean` (25,058 rxns) from
  `salvaged` (7,819 rxns) so consumers can filter to strict rows
  for mechanism-level use.
- Symmetry-group rewriting (`Build_Atom_Equivalence_Groups.py` +
  the Rebuild step's set-notation output) captures InChI equivalence
  classes for **20,887 compounds** and expresses "chemically
  indistinguishable" atoms as set notation on mapping-row endpoints,
  instead of arbitrary specific-atom picks.
- SOLR indexes the parent-level flat fields
  (`atom_mapping_data`, `atom_mapping_confidence`,
  `atom_mapping_has_symmetry_groups`) — reaction detail views and
  browse-page filters both work.

## Open considerations

### 1. Rendering symmetry sets in the UI

The `data` rows can now contain endpoints like
`cpd00012:(O#1;O#2;O#3;O#4;O#5;O#6)`. A UI that draws atom-atom
mapping lines on rendered molecular structures needs a rule for how
to visualize "atom in a set". Three patterns are viable and were
sketched in a working discussion:

- **Color-code equivalent atoms** — every atom in a group gets the
  same color; mapping is implied by color match. Breaks when a source
  set fans out to multiple destination sets (rxn00313 case: DAP's
  4 equivalent Os → Lys's 2 Os AND CO₂'s 2 Os).
- **Draw a hull around each group + hull-to-hull arrow** — works for
  the fan-out case; may visually clutter for complex reactions.
- **Interactive hover + strict-mode toggle** — recommended: default
  render uses color or hulls; hovering an atom highlights all group
  members + all possible partners; a corner toggle switches to
  "strict single-atom" mode that picks one representative per
  set-notation row.

**Action:** designate a UI-design task once the ModelSEED-UI team has
capacity. The `atom_mapping.has_symmetry_groups` boolean is the
enabling flag — pages/queries that filter it off see the same
"strict-single-atom" world as before the symmetry work.

### 2. Reactions where the missing atom is inside a broken chain

Row-level filter recovers pair rows from run-on chains but drops any
individual atom that only appeared in the cross-element portion of
the chain (typical case: a decarboxylation's carboxyl C is stitched
to the CO₂ carbon via a `C=O=O=O` chain — the chain gets split into
adjacent same-element pairs, and the C→C mapping is lost because
neither `C=O` sub-pair matches).

For rxn00313 (DAP → Lys + CO₂): the current output includes 11 rows
after symmetry rewriting, but the specific `DAP:carboxyl-C →
CO₂:C` mapping is still missing — RDT emitted it inside a chain
with element mismatches and the filter drops the cross-element parts.
Inference "the missing C is the equivalent-class member that DAP
had one of but Lys only received one from" is possible but not
implemented.

**Action:** low-priority follow-up. For biomass reachability this
doesn't matter (undirected connectivity is preserved via the O→O
rewriting). For mechanism-level ¹³C or ¹⁴C flux, users need to be
aware that C→C fates across decarboxylations may be underspecified.

### 3. Multi-molecule product stoichiometry

The `atom_mapping` row format uses compound identity (not per-molecule
instance): `cpd00012:P#1=cpd00009:P#1` means "some P of PPi maps to
some P of some Pi molecule". In rxn00001 (H₂O + PPi → 2 Pi), the
reaction produces TWO Pi molecules but the row can only reference
`cpd00009` once — the "which Pi copy" information is lost.

Symmetry rewriting compounds this: the row
`cpd00012:(P#1;P#2)=cpd00009:P#1` correctly captures that either
PPi P becomes a Pi P, but does not distinguish "the first Pi
molecule" from "the second Pi molecule".

**Action:** downstream consumers that care about per-molecule-instance
tracking must consult the reaction's stoichiometry (multiplicity
of each compound_id). No format change planned — the RDT output has
the same limitation.

### 4. RDT-flagged reactions (RDT ran but no clean rows)

1,054 reactions (57 in priority scope) had raw RDT output where every
row failed the canonical check even at the row level — element
mismatches on every row, or the whole output was chains/orphans.
These reactions have no atom_mapping ingested at all.

Sample audit hasn't been done. Some subset probably includes RDT
failures on very complex chemistry (Diels-Alder, sigmatropic
rearrangements) where the graph-alignment algorithm fundamentally
can't decide. Others are probably just bugs RDT could fix if
reproduced. Contributing back to Sebastian requires either identifying
specific reactions and asking RDT authors, or accepting the coverage
gap.

**Action:** open question. Not blocking the paper.

### 5. Pipeline-timeout reactions (~1,735)

Sebastian's rerun with an extended per-reaction budget already
recovered 354 of the original 531 priority-scope timeouts. The
remaining ~1,735 across the full DB (37 in priority scope) still
time out. Options:

- Even longer per-reaction budget (5+ min) — recovers some fraction
- Drop the `-c` (complexMode / ring exploration) flag — faster but
  loses ring-rearrangement accuracy for those reactions
- Combine both

**Action:** deferred until the local `UniversalRDT/` env has cycles
to spare (built on poplar, ready to go). Documented in the
`ATOM_MAPPING_CARBON_GAP.md`-era briefing.

### 6. Compounds without SMILES (18,621 rxns)

Fundamentally blocked on the compound-curation side — no SMILES
means RDT can't attempt the reaction. Structure-curation follow-up
(covered under Ray's `structure_update` fork and the InChIKey-mismatch
audit) may reduce this over time. Not something the mapping
pipeline can fix directly.

### 7. Wildcard SMILES (1,725 rxns)

Reactions where at least one participant has `*` in its SMILES
(unspecified R-group / abstract compound / polymer). Fundamentally
unmappable without picking a concrete placeholder atom for the
wildcard, which would fabricate structure. Sebastian committed these
to `all_rxns_with_joker.txt` for reference; excluded from clean set
by construction.

**Action:** none — permanent limitation of the "atoms as SMILES atoms"
mapping model.

### 8. Second-mapper cross-validation

Considered but not done: run a second atom-mapper (RxnMapper /
IBM RXN's transformer model, or EC-BLAST) on the same reactions and
cross-validate against RDT's output. Where they agree, confidence
compounds. Where they disagree, flag for review.

**Action:** a possible future extension. Would live in a companion
`rxnmapper_mappings.txt` in this directory + a `mapper_consensus`
field on reaction JSONs.

### 9. Symmetry-group coverage completeness

InChI's `/E:` layer only reports atom equivalences the InChI
canonicalization algorithm detected. It works for topological
symmetry (mirror planes, rotational symmetry in the Lewis structure)
but doesn't capture:

- **Tautomeric equivalence** — enol / keto interconversion swaps H
  and can make apparent atom labels swap
- **Rotamer / conformational equivalence** — freely-rotating methyl
  Hs / freely-rotating single bonds — usually the tool handles the
  common cases (methyls) but not conformational sampling
- **Charge-state equivalence** — the two Os of a carboxylate are
  equivalent by resonance; InChI does capture this for canonical
  protonation states, but not always across protonation-state changes

For the ModelSEED biochemistry these are mostly minor. If mechanistic
consumers report residual false-precision, we could layer on domain
rules (e.g., "always treat the two Os of any -COO⁻ as equivalent
regardless of what InChI says"), but this is a followup.

## Where to add new items

Append below with a `### N. Title` heading following the same shape:
one paragraph on the issue + one paragraph on **Action:**. Don't
delete closed items — mark them `**Status: RESOLVED (link to PR)**`
inline so the history stays intact.
