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
the same heuristic as the canonical reversibility, but applied to that method's
own dG:

```json
"thermodynamics": {
    "Group contribution": [4.15, 1.22, "="],
    "eQuilibrator":       [-3.46, 0.05, ">"],
    "dGPredictor":        [-3.82, 0.02, ">"]
}
```

`dGPredictor` (Wang et al. 2021; predictions staged under
`../../Biochemistry/Thermodynamics/dGPredictor/json_files/`, kJ→kcal `/4.184`)
is recorded for **every** reaction it predicts, alongside the GC/eQ records.
These per-method records sit **next to**, and never replace, the canonical
top-level `deltag` / `deltagerr` / `reversibility` fields — recording
dGPredictor does not alter the canonical free-energy value. The shared heuristic
lives in `Estimate_Reaction_Reversibility.py` (`reversibility_from_energy`); the
`Update_Reaction_*_Energies.py` scripts attach the operator when they write each
method, and `Add_Reaction_Thermodynamics_Operators.py` can (re)generate the
operators for all stored energies at any time without needing the upstream
GC / eQuilibrator / dGPredictor inputs.

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
# Backfill/refresh the per-method [energy, error, operator] triples
./Add_Reaction_Thermodynamics_Operators.py
```

These easily run together by running:
```
./Rerun_Thermodynamics.sh
```
