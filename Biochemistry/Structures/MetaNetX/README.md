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

```
$ quilt version list equilibrator/cache
0.2.0: dddb582567ce3eec2e0c984032506c7a47eabb1e7126e19cebb0bd6bce3bd9d3
0.2.1: 01e755c6234c968b73bdc061a0cfcbb2cc83ab7b6827132b2f57ef8ba6e66858
0.2.3: 6fc4e3cfb0b954488937b823c4013fe613d80644479b7a9e3fd84d524bade402
0.2.4: 43e27a997dbb7fd576b719439828d6f7402967f2530739340521574aa0b392ea
```

```
$ pip show equilibrator_api
Name: equilibrator-api
Version: 0.2.5
Summary: Standard reaction Gibbs energy estimation for biochemical reactions.
Home-page: https://gitlab.com/elad.noor/equilibrator-api
Author: Elad Noor
Author-email: noor@imsb.biol.ethz.ch
License: MIT License
Location: /anaconda2/envs/msd_env/lib/python3.7/site-packages
Requires: numpy, pyparsing, sbtab, component-contribution, matplotlib, pandas, nltk, scipy, equilibrator-cache, optlang
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

```
Download MetaNetX file:

https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv
(Downloaded version 2019/02/13)

```


## Latest instructions

### Clear everything out
```
pip uninstall equilibrator_api equilibrator_cache component_contribution
quilt rm equilibrator/cache
```

### Re-install latest code
```
pip install --pre equilibrator_api
```

### To find SQLite table, need to edit code to print out new location of table

```
Edit /anaconda2/envs/py3/lib/python3.6/site-packages/equilibrator_cache/api.py to print "location"
python
>>> from equilibrator_api import ccache
$ sqlite3 /var/folders/fx/lnfg_9152tdg2ffl8236g0gr0000gq/T/tmp5sffp_kd/compounds.sqlite
```

```
sqlite> .headers on
sqlite> .mode tabs
sqlite> .output cpds.dump
sqlite> select ci.accession,c.inchi_key,c.smiles from compounds as c,compound_identifiers as ci where c.id==ci.compound_id and ci.registry_id==4 order by c.id,ci.accession;
sqlite> .quit
```