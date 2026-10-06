# Agent instructions for ModelSEEDDatabase

## Never read large JSON/TSV data files directly into context

This repository's `Biochemistry/` tree holds the live database: ~60 reaction
shards and ~50 compound shards (each several MB of JSON), plus large flat
files like `Unique_ModelSEED_Structures.txt` (13+ MB) and grading tables like
`reaction_grades.tsv` (tens of thousands of rows). Do **not** use the `Read`
tool (or `cat`/`head -c large`) to load one of these files wholesale into the
conversation — it burns an enormous amount of context for data a script can
answer in one line, and large files get silently truncated by the Read tool
regardless.

**Instead:** write a short, throwaway Python (or `jq`/`awk`) snippet via the
Bash tool, run it, and only bring the *result* (a count, a table, a handful of
rows) into context. Treat every file under `Biochemistry/`,
`Biochemistry/Thermodynamics/SourceGrading/results/`, and similar data
directories this way — load it in the script, query/aggregate it there, print
only the answer.

Examples:

```bash
# WRONG: reads the whole 5MB shard into the conversation
# (Read tool on Biochemistry/reaction_00.json)

# RIGHT: load it in a script, print only what's needed
python3 -c "
import json
d = json.load(open('Biochemistry/reaction_00.json'))
print(len(d), 'reactions')
print(d[0]['id'], d[0].get('thermo-evidence'))
"
```

```bash
# Aggregating across all reaction shards without ever printing a record
python3 -c "
import json, glob
from collections import Counter
c = Counter()
for f in sorted(glob.glob('Biochemistry/reaction_*.json')):
    for r in json.load(open(f)):
        e = r.get('thermo-evidence')
        if e:
            c[e['grade']] += 1
print(c)
"
```

```bash
# TSV: use csv.DictReader in a script rather than grepping/catting the file
python3 -c "
import csv
from collections import Counter
c = Counter()
for r in csv.DictReader(open('Biochemistry/Thermodynamics/SourceGrading/results/thermo_grades/reaction_grades.tsv'), delimiter='\t'):
    c[r['best_grade']] += 1
print(c)
"
```

If you need to inspect a handful of specific records, filter for them inside
the script (`if r['id'] in (...)`) and print just those — never dump the
whole structure to find them.

This applies to `grep`/`ripgrep` too when the match count could be large:
prefer `grep -c` or piping through a script to count/summarize rather than
printing every matching line of a multi-thousand-line file.
