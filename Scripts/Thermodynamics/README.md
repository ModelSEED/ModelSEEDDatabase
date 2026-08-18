# ModelSEED Biochemistry Database Scripts (Thermodynamics)

Here we have the scripts used to handle the thermodynamics data that
came from both the Group Contribution approach (<a
href="https://doi.org/10.1529/biophysj.107.124784">Jankowski et
al. 2008</a>) and from eQuilibrator (<a
href="https://doi.org/10.1371/journal.pcbi.1003098">Noor et
al. 2013</a>). Our approach in handling this data is described in the
<a
href="https://www.biorxiv.org/content/10.1101/2020.03.31.018663v2">paper</a>.

## Order of execution

The general order is that the energies from the application of the Group Contribution (GC) approach
are stored in the database first, and then the energies from eQuilibrator (EQ), which, in most
cases, take precedence, are used to overwrite the energies in the database

### Additive per-method thermodynamics

Each reaction keeps every method's estimate **additively** in its
`thermodynamics` dict rather than collapsing them into a single value. Each
method holds an `[energy, error, operator]` triple, where the operator (`>`,
`<`, `=`, or `?`) is that estimate's own thermodynamic direction — computed with
that method's own rule set (see *Source-specific reversibility heuristics*
below) applied to that method's own dG:

```json
"thermodynamics": {
    "Group contribution":    [4.15, 1.22, "="],
    "eQuilibrator":          [-3.46, 0.05, ">"],
    "dGPredictor":           [-3.82, 0.02, ">"],
    "dGPredictor-ModelSEED": [-3.77, 0.87, ">"]
}
```

`dGPredictor` (Wang et al. 2021; predictions staged under
`../../Biochemistry/Thermodynamics/dGPredictor/json_files/`, kJ→kcal `/4.184`)
is recorded for **every** reaction it predicts, alongside the GC/eQ records.
`dGPredictor-ModelSEED` is the same dGPredictor model **retrained on the
ModelSEED compound structures** (expanded group vocabulary; staged under
`../../Biochemistry/Thermodynamics/dGPredictor/modelseed_retrained_dG.json`,
kJ→kcal `/4.184`), recorded as its **own** additive method for the 31,924
reactions it predicts — next to, and never replacing, the original KEGG-based
`dGPredictor` record. These per-method records sit **next to**, and never
replace, the canonical top-level `deltag` / `deltagerr` / `reversibility`
fields — recording dGPredictor does not alter the canonical free-energy value. The rule sets
live in `reversibility_heuristics.py`, reached via
`Estimate_Reaction_Reversibility.reversibility_from_energy(..., source=<label>)`; the
`Update_Reaction_*_Energies.py` scripts attach the operator when they write each
method, and `Add_Reaction_Thermodynamics_Operators.py` can (re)generate the
operators for all stored energies at any time without needing the upstream
GC / eQuilibrator / dGPredictor inputs.

### Source-specific reversibility heuristics

The direction cascade lives in `reversibility_heuristics.py` and is selected
**per thermodynamic data source**, because the sources do not fail the same way.
`Estimate_Reaction_Reversibility.py` picks the rule set from its argument:

```
./Estimate_Reaction_Reversibility.py            # top-level deltag, GC rules
./Estimate_Reaction_Reversibility.py GC         # Group contribution, GC rules
./Estimate_Reaction_Reversibility.py EQ         # eQuilibrator energies, EQ rules
./Estimate_Reaction_Reversibility.py EQ --heuristics EQ2   # eQuilibrator 2.0 rules
./Estimate_Reaction_Reversibility.py EQ --heuristics GC    # pre-split behaviour
```

**`GC` is the default.** Any source without a rule set of its own — including
both dGPredictor variants — falls back to it, and the GC cascade is unchanged,
so `GC` and unfiltered runs still reproduce the historical report exactly.

| set | cascade | citation |
|-----|---------|----------|
| `GC` | ATPS → ABC transporter → MdeltaG bounds → mMdeltaG band → low-energy points → `=` | Jankowski et al. 2008 |
| `EQ` | ATPS → ABC transporter → undecomposable → uncorrected transport → ln Γ ± 1σ → `=` | Beber et al. 2022 over Noor et al. 2012 |
| `EQ2` | same, with the confidence margin off (ln Γ point estimate) | Flamholz et al. 2012 + Noor et al. 2012 |

The same registry drives the per-method operators, so
`thermodynamics['eQuilibrator'][2]` is now computed with the EQ rules while the
other methods keep the GC rules.

#### What the EQ rules do

Direction comes from the **reversibility index** Γ of Noor et al. 2012 — the
fold change every reactant concentration must undergo to reverse the reaction:

```
ln Γ = (2 / Σ|ν|) · ΔG′m / RT        ΔG′m = ΔG′° + RT · Σν · ln(10⁻³)
```

with water and protons excluded from both sums, matching
`equilibrator_cache.reaction.items(protons=False, water=False)`. A reaction is
irreversible when `|ln Γ| − σ > ln(1000)`, i.e. Noor's headline 3 µM–3 mM window
around 100 µM. `Scripts/Tests/test_eq_heuristics.py` verifies this reproduces
eQuilibrator's own published `ln_reversibility_index` (column 4 of
`MetaNetX_Reaction_Energies.tbl`) on 17,771 reactions.

Beber 2022 defines no reversibility index; what eQuilibrator 3.0 contributes is
rigorous uncertainty, and that is what the two gates ahead of the index rule
encode. Both address defects in the stored eQuilibrator data:

* **`EQ:undecomposable`** — a σ of ~1e5 kJ/mol (`RMSE_inf`) marks a reaction
  eQuilibrator could not decompose. For the majority of these it does not return
  a prediction at all: `GibbsEnergyPredictor.standard_dg` short-circuits with

  ```python
  if residual:
      return Q_(0, "kJ/mol").plus_minus(self.preprocess.RMSE_inf)
  ```

  i.e. **literally zero**, and the nonzero ΔG′° we store is only the Legendre/pH
  transform that `standard_dg_prime` adds on top of that zero. Of the 4,607
  σ-flagged reactions in the table, **2,895** have σ exactly `RMSE_inf` and are
  this zeroed case — 61 of them are still exactly `0.000000` because they have no
  net proton change for the transform to act on. The remaining **1,712** carry
  `RMSE_inf * ||sigma_inf||` (e.g. 1/√2) and do have a real mean, but with an
  undetermined direction in the covariance.

  Either way the value cannot support a directional call, so the gate returns
  `?`. Previously the GC bounds rule could never fire at that width and these
  fell through to a permissive `=`. (Real σ tops out at 65.35 kcal/mol against a
  marker of 23,900.57, so the cut at 1e4 kJ/mol sits in an empty gap.)
* **`EQ:transport-uncorrected`** — `Retrieve_eQuilibrator_Reactions_Energies.py`
  keys its MetaNetX formula on compound id and therefore discards compartment,
  so any species on both sides nets out; 1,102 transport reactions carry a ΔG′°
  for a *different* reaction. Independently, Beber 2022 notes the transformed
  framework needs a `−N_H·RT·ln(10^ΔpH) − Q·FΔΦ` term across a membrane that
  this pipeline never applies. ATP synthase and ABC transporters are still
  decided structurally, ahead of this gate.

Reactions whose ln Γ interval merely straddles the threshold are reported as
`EQ:ambiguous` and called `=`: Noor treats Γ as a continuous index and reserves
the directional call for clear cases.

#### Energy source for the EQ run

`Estimate_Reaction_Reversibility.py EQ` reads
`thermodynamics['eQuilibrator']`'s own dG and σ, **not** the top-level `deltag`.

Reading `deltag` was correct as originally designed. From 2019 (`c263e233`)
through 2023-09-11 (`17c9739b`), `Update_Reaction_eQuilibrator_Energies.py`
overwrote the canonical energy and stamped the note that gates this step:

```python
reactions_dict[rxn]['deltag']    = float(eq_reactions[rxn]['dg'])
reactions_dict[rxn]['deltagerr'] = float(eq_reactions[rxn]['dge'])
if('EQU' not in notes_list): notes_list.append('EQU')
```

so `deltag` *was* the eQuilibrator value wherever `EQU` was set, and
`DB_LEVEL_NOTE = {"EQ": "EQU"}` gated it exactly right.

**`3e50646b` (2023-09-13) removed the write half and left the read half in
place.** The updater moved to `thermodynamics['eQuilibrator']` only, stopping
both the `deltag` overwrite and the `EQU` stamp — but this step kept reading
`deltag`. Later GC rebuilds (Convention A, `ad34d6ab`) and the PR #265
promotion then rewrote `deltag` from other sources. Today only **1,797 of
25,028** reactions with an eQuilibrator record have
`deltag == thermodynamics['eQuilibrator'][0]`; restricted to the 17,094 still
carrying `EQU`, only **1,586** do. That note is now a fossil of the pre-2023
pipeline. (`EQU` is not `EQC`: the retrieval script writes `EQC`/`EQP` to mean
"all reagents had eQuilibrator structures", a different claim.)

Pointing the EQ run at the eQuilibrator sublist restores the original intent.
The EQ rules also need eQuilibrator's own σ, which the top-level `deltagerr`
never carries.

Consequence: canonical `reversibility` for these reactions now derives from the
eQuilibrator energy while canonical `deltag` mostly remains the GC value.
Reconciling `deltag` itself is the job of
`Promote_Reaction_Thermodynamics_to_Canonical.py`, which is deliberately not
part of `Rerun_Thermodynamics.sh`.

The underlying thermodynamics data is kept in
`../../Biochemistry/Thermodynamics`. The decomposition of molecular
structures and their resulting energies for both the older group
contribution approach and the newer eQuilibrator approach are stored in
the `ModelSEED` and `eQuilibrator` directories.

As an addendum, the two scripts used to update the energies from
eQuilibrator are in this folder, but are dependent on files in
`../../Biochemistry/Structures/MetaNetX`:
```
./Retrieve_eQuilibrator_Compound_Energies.py
./Retrieve_eQuilibrator_Reactions_Energies.py
```

If the underlying thermodynamics data in `../../Biochemistry/Thermodyanmics` hasn't changed,
then running these six commands should not cause any changes to appear in the database.

```
./Update_Compound_GroupContribution_Energies.py
./Update_Reaction_GroupContribution_Energies.py
./Estimate_Reaction_Reversibility.py GC
./Update_Compound_eQuilibrator_Energies.py
./Update_Reaction_eQuilibrator_Energies.py
./Estimate_Reaction_Reversibility.py EQ
# Record dGPredictor additively for every predicted reaction (no canonical change)
./Update_Reaction_dGPredictor_Energies.py
# Record the ModelSEED-retrained dGPredictor as its own additive method
./Update_Reaction_dGPredictor_ModelSEED_Energies.py
# Backfill/refresh the per-method [energy, error, operator] triples
./Add_Reaction_Thermodynamics_Operators.py
```

These easily run together by running:
```
./Rerun_Thermodynamics.sh
```
