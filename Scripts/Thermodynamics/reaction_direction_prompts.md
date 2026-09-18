# Reaction-direction prompts for LLM classification

Two prompts for asking a language model to classify the *in vivo* direction of a
ModelSEED reaction.

**Prompt A** is the one already in production, reproduced verbatim. **Prompt B**
is a proposed extension that additionally supplies the thermodynamic estimates
this database computes, so the model is asked to reconcile biology with data
rather than to work from the equation alone.

Prompt A is deliberately unchanged so that A-vs-B is a clean comparison: the
only difference is the presence of the thermodynamic block.

---

## Provenance of Prompt A

`AICurationUtils.analyze_reaction_directionality`, in
`kbutillib/ai_curation_utils.py`. Driven by `regen_directionality.py`, answered
by Claude Opus 4.8 through the Argonne Argo gateway, and cached to
`AICurationCacheReactionDirectionality.json`. Downstream consumers are
`build_ai_direction_variant.py` (turns the labels into a reversibility variant)
and `estimate_directions_literature.py` (compares them against Jankowski 2008,
Flamholz 2012 and Beber 2022).

Three behaviours of the harness matter when interpreting results, and are not
visible in the prompt text:

1. **Reactions may be inverted before asking.** Where the harness reverses a
   reaction for presentation, it flips `forward` / `reverse` back afterwards and
   appends a note to `other_comments`. Cached directions are post-correction,
   not raw model output.
2. **The cache is never invalidated.** A reaction already present is not
   re-queried, so a cache file may mix results from different model versions or
   prompt revisions unless it was regenerated wholesale.
3. **Reactions whose ID prefix is a utility prefix are skipped** and return
   `None`.

---

## Prompt A — as used (neutral, no data)

### System

```text
You are an expert in biochemistry and molecular biology. 
You will receive a biochemical reaction and must evaluate it for stoichiometric 
correctness and biological directionality.

Respond strictly in valid JSON with **no text outside the JSON**. 
All keys and string values must use double quotes. 
Use only plain ASCII characters.
```

### User

```text
Analyze the following reaction for stoichiometric correctness and 
directionality in vivo. 

Return a JSON object in this exact format:

{
"errors": ["error 1", "error 2"],
"directionality": "forward|reverse|reversible|uncertain",
"other_comments": "Brief general comments about the reaction so I know you understood the input.",
"confidence": "high|medium|low|none"
}

Reaction:
<reaction string>
```

---

## Prompt B — proposed, with thermodynamic estimates

Same task and same output contract, plus the per-source estimates. The added
fields are additive so that A and B remain directly comparable on
`directionality` and `confidence`.

### System

```text
You are an expert in biochemistry and molecular biology.
You will receive a biochemical reaction, together with independent computational
estimates of its standard transformed Gibbs free energy of reaction, and must
evaluate the reaction for stoichiometric correctness and biological
directionality.

The estimates are predictions, not measurements. They disagree with each other,
they carry different and imperfectly calibrated uncertainties, and none of them
knows anything about the physiological concentrations that actually set the
direction in a cell. Treat them as evidence to be weighed against biochemical
knowledge, not as an answer to be reported back.

Respond strictly in valid JSON with **no text outside the JSON**.
All keys and string values must use double quotes.
Use only plain ASCII characters.
```

### User

```text
Analyze the following reaction for stoichiometric correctness and
directionality in vivo.

You are given up to three independent estimates of the standard transformed
Gibbs free energy of reaction, dGr'0, in kcal/mol at pH 7.0, ionic strength
0.25 M, pMg 3.0 and 298.15 K. Each is reported with its own uncertainty.
A value may be absent if that method could not estimate this reaction.

  Group contribution   a group-decomposition method; broad coverage, and its
                       stated uncertainty is known to be conservative
  dGPredictor          a machine-learning model over molecular fragments
  eQuilibrator         component contribution, fitted to experimental data;
                       its stated uncertainty is known to be optimistic and
                       excludes any contribution from the protonation model

Guidance on how to use them:

- A reaction is conventionally treated as effectively irreversible when
  |dGr'0| is greater than about 2 kcal/mol, and as reversible within that band.
  This is a convention, not a law: physiological concentrations can drive a
  thermodynamically unfavourable reaction forward, and mass action routinely
  overrides a modest dGr'0.
- Where the estimates agree, they may still be wrong together; two of these
  methods share group-contribution lineage.
- Where they disagree, say so and say which you find more credible and why.
- If the biochemistry contradicts all three estimates, prefer the biochemistry
  and state that you are doing so.

Return a JSON object in this exact format:

{
"errors": ["error 1", "error 2"],
"directionality": "forward|reverse|reversible|uncertain",
"other_comments": "Brief general comments about the reaction so I know you understood the input.",
"confidence": "high|medium|low|none",
"thermodynamic_assessment": "Brief note on whether the estimates support your answer, and how you resolved any disagreement between them.",
"agrees_with_thermodynamics": "yes|no|partly|no data"
}

Reaction:
<reaction string>

Thermodynamic estimates (dGr'0, kcal/mol):
  Group contribution : <value> +/- <error>
  dGPredictor        : <value> +/- <error>
  eQuilibrator       : <value> +/- <error>
```

### Worked example

`rxn00001`, diphosphate phosphohydrolase, where the sources disagree in sign:

```text
Reaction:
(1) H2O[0] + (1) PPi[0] => (2) Phosphate[0] + (1) H+[0]

Thermodynamic estimates (dGr'0, kcal/mol):
  Group contribution : 4.18 +/- 2.24
  dGPredictor        : -3.78 +/- 1.48
  eQuilibrator       : -4.06 +/- 0.05
```

A second useful case is `rxn00549`, fructose-1,6-bisphosphatase, where Group
contribution says +7.22 ± 3.82 and eQuilibrator says −3.00 ± 0.17 — a reaction
whose physiological irreversibility is not in doubt, so it tests whether the
model defers to biology when a source contradicts it.

---

## Notes for whoever runs this

**Units.** ModelSEED stores these in kcal/mol. eQuilibrator's own API reports
kJ/mol; do not mix them.

**Sentinels must be stripped before templating.** Group contribution stores
`10000000.0` for "no estimate" on 26,555 of 56,012 reactions, and eQuilibrator
uses very large uncertainties for the same purpose. Passing either through as a
number would be worse than omitting the source.

**Direction vocabulary differs from the database.** The prompt returns
`forward|reverse|reversible|uncertain`; ModelSEED stores `>`, `<`, `=`, `?`. The
mapping happens downstream in `build_ai_direction_variant.py`.

**The uncertainty characterisations in Prompt B are measured, not folklore.**
Against 802 experimentally-anchored reactions, Group contribution overstates its
error by roughly 2.2x, dGPredictor by 1.5x, and eQuilibrator understates by
1.6x. Stating this in the prompt is a deliberate choice: it may improve
weighting, but it also steers the model, so an A/B comparison should ideally
include a variant with that guidance removed.

**Suggested comparison.** Run A and B over the same reaction set and compare
`directionality` agreement, agreement with the reversibility rule applied to
each source's dGr'0, and agreement with the gold-tier subset of
`Biochemistry/Thermodynamics/SourceGrading/results/thermo_grades/`, where an
experimental measurement exists.
