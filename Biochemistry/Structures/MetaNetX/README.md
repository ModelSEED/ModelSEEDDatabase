## Loading eQuilibrator cache

The "latest" set of data in the eQuilibrator cache does not work with equilibrator_api version 0.2.0
So, here is a means to make sure that you retrieve the cache from quilt that will work:

```
bash-3.2$ pip install quilt
bash-3.2$ quilt install equilibrator/cache -x 7a100a
```

## Dumping structures from eQuilibrator cache

In `./Find_MNX_Structures_in_MS.pl` We cross-check the available mappings between
MetaNetX identifiers and InChIKeys to make sure that we use only the mappings that still
exist in MetaNetX, but to do so, we dump the structures from the eQuilibrator cache.

You can use the eQuilibrator API directly to recall the structures, but it takes
more than 3 hours on a laptop, and this process takes a few minutes.

### Step 1

The first time you use the equilibrator_api, it will populate an sqlite table with the data in quilt,
so here, we force that to happen

```
bash-3.2$ pip install equilibrator_api
bash-3.2$ python
>>> from equilibrator_api import ccache
```

### Step 2

Secondly, you have to find where the sqlite table is on your system. Here we find it within our
conda environment, but you may have to use a different path depending on how your python is configured
```
bash-3.2$ find /anaconda2/envs -name compounds.sqlite
/anaconda2/envs/py3/lib/python3.6/site-packages/equilibrator_cache/cache/compounds.sqlite
```

### Step 3

Finally, you load the table in sqlite, and dump:

```
bash-3.2$ sqlite3 /anaconda2/envs/py3/lib/python3.6/site-packages/equilibrator_cache/cache/compounds.sqlite
SQLite version 3.27.2 2019-02-25 16:06:06
Enter ".help" for usage hints.
sqlite> .headers on
sqlite> .mode tabs
sqlite> .output eq_cpds.dump
sqlite> select mnx_id,inchi_key from compounds;
sqlite> .quit
```
