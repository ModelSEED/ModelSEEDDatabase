## Order of execution

The general order is that the energies from the application of the Group Formation (GF) approach
are stored in the database first, and then the energies from eQuilibrator (EQ), which, in most
cases, take precedence, are used to overwrite the energies in the database

The underlying thermodynamics data is kept in
`../../Biochemistry/Thermodynamics`. The decomposition of molecular
structures and their resulting energies for both the older group
formation approach and the newer eQuilibrator approach are stored in
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
./Update_Compound_GroupFormation_Energies.py
./Update_Reaction_GroupFormation_Energies.py
./Estimate_Reaction_Reversibility.py GF
./Update_Compound_eQuilibrator_Energies.py
./Update_Reaction_eQuilibrator_Energies.py
./Estimate_Reaction_Reversibility.py EQ
```
