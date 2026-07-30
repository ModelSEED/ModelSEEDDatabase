# Protein-carrier cofactor census — 2026-07-29

Census of candidate compound classes for §5's database-wide protein-carrier cofactor standardisation, run against `dev` HEAD **`bbbedbe`**.

## Method

Compounds are pattern-matched against a per-class list of name substrings (case-insensitive, matching both the primary `name` field and the `Name:` alias line in `aliases`). Results are cross-tabbed against the existing `acps_formula_charge.tsv` override set to distinguish already-covered from open work.

**Caveats:** name-pattern matching is over-broad. Some hits are compounds with the substring but no carrier semantics (e.g. cpd09506 "(S)-ACPA" is a specific small molecule, not an acyl-carrier form). A follow-up pass must apply a structural filter (compounds with `R` in the formula, or with a `*` in the SMILES that represents a protein attachment) before overrides are proposed.

## Class-by-class

### Acyl-ACP family

- **Name patterns:** `-[acp]`, `[acp]`, `-acp`, `acyl-carrier`, `acyl carrier`
- **Total candidates:** 830 compounds
- **Already have a pantetheine-inclusive override** (row in `acps_formula_charge.tsv`): **629** (76%)
- **No override yet:** **260** compounds. Some of these are true acyl-ACPs missed by the current curation pass; some are false positives (name substring, no carrier semantics).

Next round of §5 work on ACP: apply the structural filter to the 260 open candidates, run `Compute_ACP_Overrides.py` on the survivors to propose formula rows, apply `Report_Formula_Change_Impact.py` to filter down to the subset that would not newly imbalance any reaction (identical procedure to the earlier bulk-add-47-safe-overrides commit), commit as a bulk-add PR.

### Biotinyl-carrier family

- **Name patterns:** `biotinyl-`, `-biotinyl`, `bccp`, `biotin-carboxyl`, `carboxybiotin`, `carboxy-biotin`, `holo-bccp`, `apo-bccp`
- **Total candidates:** 21 compounds

Manual inspection of the top hits:

| cpd_id | name | formula | Carrier? |
|---|---|---|---|
| cpd12543 | Carboxybiotin-carboxyl-carrier protein | C18H26N5O6R2S/-1 | Yes (has R groups) |
| cpd12848 | Holo-[carboxylase] | C17H27N5O4R2S/0 | Yes (has R groups) |
| cpd01302 | Biotinyl-CoA | C31H46N9O18P3S2/-4 | No (CoA-bound, not protein) |
| cpd03293 | Biocytin | C16H28N4O4S/0 | Border — biotin-lysine adduct, free small molecule |
| cpd03517 | Biotinyl-5-AMP | C20H27N7O9PS/-1 | No (activated form, not carrier) |
| cpd11415 | carboxybiotin | C11H14N2O5S/-2 | No (free small molecule) |

**True protein-carrier subset is small — approximately 4-5 compounds.** Small enough to hand-curate directly.

Standardisation proposal for the biotinyl class: for every carrier-form biotinyl compound (has `R` in formula, protein attachment via lysine ε-amine amide bond), include the biotin moiety in the stored formula. Because the class is small, a single hand-authored PR extending the override framework covers it. (An override file for biotinyl carriers may be `Biochemistry/Curation/overrides/biotinyl_formula_charge.tsv`, mirroring the ACP file's schema.)

### Lipoyl-carrier family

- **Name patterns:** `lipoyl-`, `-lipoyl`, `lipoamide`, `lipoylprotein`, `dihydrolipoyl`, `lipoyl protein`, `lipoyl-h-protein`, `lipoyl-e2`
- **Total candidates:** 82 compounds

The vast majority of these are named as free lipoamide or S-acyl-dihydrolipoamide compounds (e.g. cpd00213 Lipoamide, cpd00449 Dihydrolipoamide, cpd00836 S-Acetyldihydrolipoamide). **None of the first 6 examples have `R` in the formula.** This is the boundary problem the working session identified: the DB stores these as small-molecule lipoamide forms, but the enzymatic reactions (PDH, KGDH, BCKDH, glycine cleavage H-protein) genuinely represent lipoic acid *covalently bound to a lysine ε-amine on a protein*. The current representation therefore either (a) mass-balances against a fictional free-lipoamide chemistry, or (b) does not participate in reactions that require protein-bound loading.

**Standardisation for the lipoyl class is substantially more work than for biotinyl or ACP** because it may require redrawing SMILES structures to add the protein-side `R` group, not just changing the formula field. A first pass should audit the 82 compounds to determine:

1. Which are genuine free small molecules and should remain unchanged (lipoic acid, dihydrolipoic acid).
2. Which represent protein-bound forms and should be restructured to include an `R` protein attachment plus the corresponding formula override.

## Candidate extensions (deferred per plan)

- **Covalent FMN/FAD** (e.g. histidyl-FAD of succinate dehydrogenase): need a census pass.
- **Molybdopterin cofactors**: need a census pass.
- **Heme c** (thioether-attached to cytochromes): need a census pass.

## Summary counts

| Class | Candidates | Already overridden | Open | Standardisation effort |
|---|---:|---:|---:|---|
| Acyl-ACP | 830 | 629 | 260 (with structural filter, actual might be smaller) | Extend existing framework via bulk-add PR |
| Biotinyl carrier | 21 | 0 | ~4-5 true carriers after inspection | Small hand-authored PR |
| Lipoyl carrier | 82 | 0 | Requires per-compound audit | May need SMILES restructuring, not just formula override |

## Decisions to close for §5

1. Which of the 260 open ACP candidates are true acyl-ACPs after the structural filter? Run the filter, propose overrides via `Compute_ACP_Overrides.py`, gate through `Report_Formula_Change_Impact.py`.
2. Add a biotinyl override file (or extend `acps_formula_charge.tsv` to hold both classes)?
3. Scope the lipoyl audit: full 82 compounds, or start with the compounds actually referenced by priority-scope reactions?
4. Defer covalent FMN/FAD, molybdopterin, heme c to a future paper unless a quick census shows they touch a material number of priority reactions.
