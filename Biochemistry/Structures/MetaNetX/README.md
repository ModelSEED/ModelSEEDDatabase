## Dumping structures from eQuilibrator cache

In `./Find_eQuilibrator_Structures_in_ModelSEED.py` We cross-check the available mappings between
InChIKeys in eQuilibrator and ModelSEED, but to do so, we dump the structures from the eQuilibrator cache.

NB: previously we also checked with the latest version of MetaNetX, but the versioning of the MetaNetX
identifiers has made this untenable, so we check directly with eQuilibrator.

### Step 1

We install eQuilibrator's cache directly using `quilt`:
```
% quilt install equilibrator/cache
```

At time of writing, we used this version of the data from equilibrator:
```
$ quilt version list equilibrator/cache | tail -n1
0.2.10: d11d9ebd1cf25d95d8c1369a8ca2bbe46eb734f8b2b9a4225f16de614104e12c
```

### Step 2

We export the cached compound structures:
```
% quilt export equilibrator/cache equilibrator
```

### Step 3

The exported structures are in SQLite format so we load them and re-dump
them in a flat-file format:
```
% sqlite3 equilibrator/compounds.sqlite
SQLite version 3.39.3 2022-09-05 11:02:23
Enter ".help" for usage hints.
sqlite> .headers on
sqlite> .mode tabs
sqlite> .output eq_cpds.tsv
sqlite> select ci.accession,c.inchi_key,c.smiles from compounds as c,compound_identifiers as ci where c.id==ci.compound_id and ci.registry_id==4 order by c.id,ci.accession;
sqlite> .quit
```

### Step 4

We run the script to compare the two set of structures:
```
% ./Find_eQuilibrator_Structures_in_ModelSEED.py
```

You can then compare how the set of available structures has changed:
```
git diff Structures_in_ModelSEED_and_eQuilibrator.txt
```
