## Order of execution

The general order is that the energies from the application of the Group Formation (GF) approach
are stored in the database first, and then the energies from eQuilibrator (EQ), which, in most
cases, take precedence, are used to overwrite the energies in the database

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

As an addendum, since we keep these in the repository, these are the two scripts used to update
the energies from eQuilibrator:
```
./Retrieve_eQuilibrator_Compound_Energies.py
./Retrieve_eQuilibrator_Reactions_Energies.py
```