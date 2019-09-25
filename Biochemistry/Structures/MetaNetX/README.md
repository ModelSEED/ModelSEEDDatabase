## Dumping structures from eQuilibrator cache

In `./Find_MNX_Structures_in_MS.pl` We cross-check the available mappings between
MetaNetX identifiers and InChIKeys to make sure that we use only the mappings that still
exist in MetaNetX, but to do so, we dump the structures from the eQuilibrator cache.

### Step 1a

The first time you use the equilibrator_api, it will populate an sqlite table with the data in quilt,
so here, we force that to happen.

```
bash-3.2$ pip install equilibrator_api
bash-3.2$ python
>>> from equilibrator_api import ccache
```

At time of writing, we used this version of the data from equilibrator:
```
$ quilt version list equilibrator/cache
0.2.0: dddb582567ce3eec2e0c984032506c7a47eabb1e7126e19cebb0bd6bce3bd9d3
0.2.1: 01e755c6234c968b73bdc061a0cfcbb2cc83ab7b6827132b2f57ef8ba6e66858
0.2.3: 6fc4e3cfb0b954488937b823c4013fe613d80644479b7a9e3fd84d524bade402
0.2.4: 43e27a997dbb7fd576b719439828d6f7402967f2530739340521574aa0b392ea
```

along with this version of equilibrator_api:
```
$ pip show equilibrator_api
Name: equilibrator-api
Version: 0.2.5
```

### Step 1b

The `./Find_MNX_Structures_in_MS.pl` also uses a file easily downloaded from MetaNetX:
`https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv`

At time of writing, we downloaded MetaNetX version: `2019/02/13`

### Step 2a

To load the sqlite table, you first have to find where the sqlite table is on your system. Here we find it within our
conda environment, but you may have to use a different path depending on how your python is configured and version of
equilibrator_api
```
bash-3.2$ find /anaconda2/envs/msd-env -name compounds.sqlite
/anaconda2/envs/msd-env/lib/python3.6/site-packages/equilibrator_cache/cache/compounds.sqlite
```

### Step 2b

You then load the table in sqlite, and dump the relevant fields:

```
bash-3.2$ sqlite3 /anaconda2/envs/msd-env/lib/python3.6/site-packages/equilibrator_cache/cache/compounds.sqlite
SQLite version 3.27.2 2019-02-25 16:06:06
Enter ".help" for usage hints.
sqlite> .headers on
sqlite> .mode tabs
sqlite> .output eq_cpds.dump
sqlite> select ci.accession,c.inchi_key,c.smiles from compounds as c,compound_identifiers as ci \
	where c.id==ci.compound_id and ci.registry_id==4 order by c.id,ci.accession;
sqlite> .quit
```

### Step 3

Finally, you can run `./Find_MNX_Structures_in_MS.pl` and use git to check whether any of the structures have been updated:
```
git diff Structures_in_ModelSEED_and_eQuilibrator.txt
```
