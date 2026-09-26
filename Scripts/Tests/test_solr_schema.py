#!/usr/bin/env python3
"""Guard that the SOLR schemas and the compiler still agree.

Neither core declares a `dynamicField`, so the schemas are strict: Solr rejects
a document carrying a field the schema does not know. That makes the compiler
and the schema a matched pair, and nothing enforced it. They drifted -- the
multi-source pKa work and the evidence-grading work both taught
`Compile_Biochemistry_for_SOLR.py` to emit new fields (pka_kind, pka_number,
pkb_number, grade, assessment, cross_source, evidence_*) without declaring any
of them, which would have failed on every affected document at index time:
67,933 pKa children and 33,099 reactions.

Six classes of drift are checked:

1. **Malformed schema** -- the XML no longer parses.
2. **Duplicate declarations** -- the same field declared twice, where the
   second silently wins.
3. **Undeclared fields** -- the compiler emits something the schema lacks.
   This is the failure that motivated the file.
4. **Type disagreement** -- a value that will not survive the declared type,
   e.g. a string in a `float` field.
5. **Cardinality disagreement** -- a list emitted into a single-valued field.
6. **Key integrity** -- the uniqueKey is present, non-empty and unique across
   parents and children, and every child carries a `doc_type`.
7. **Refusal sentinels** -- REPORTED, NOT FAILED. Group contribution writes
   dg = 1e7 and eQuilibrator writes sigma >= 2500 kcal/mol; both mean "no
   estimate". Indexing them was a deliberate decision (2026-09-10), so the
   count is printed as a watch figure rather than an error. It is worth
   watching because a sort or facet on energy will put these at the extreme,
   and the paper excludes exactly these records from every count -- so index
   totals and manuscript totals will not agree by construction.

The nested-document container labels (`pkas`, `thermodynamics`,
`stoichiometry`, `thermo_evidence`) are deliberately NOT declared as fields:
Solr routes labelled child documents through `_nest_path_`, which the schemas
do declare. They are whitelisted here, and the whitelist is itself checked --
a new container name has to be added consciously rather than slipping through.

Usage:  python3 test_solr_schema.py [--compile]
        --compile regenerates the solr_*.json first (slow, ~1 min); without it
        the existing compiler output is used and its age is reported.
"""
import json
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SOLR = os.path.join(ROOT, 'Solr')
CORES = {
    'compounds': os.path.join(SOLR, 'compilation', 'solr_compounds.json'),
    'reactions': os.path.join(SOLR, 'compilation', 'solr_reactions.json'),
}
# Child-document container labels, routed by _nest_path_ rather than declared.
NEST_CONTAINERS = {
    'compounds': {'pkas', 'thermodynamics'},
    'reactions': {'stoichiometry', 'thermodynamics', 'thermo_evidence'},
}
NUMERIC = {'int': int, 'long': int, 'float': float}
# Refusal markers, not energies. Mirrors Scripts/Thermodynamics: GC writes
# SENTINEL_DG = 1e7, eQuilibrator writes sigma >= 2500 kcal/mol.
SENTINEL_DG = 1.0e7
SENTINEL_SIGMA = 2500.0


def declared_fields(core):
    """{name: {'type':..., 'multi':bool}} plus the raw list, for duplicates."""
    path = os.path.join(SOLR, 'cores', core, 'schema.xml')
    root = ET.parse(path).getroot()
    names, out = [], {}
    for f in root.iter('field'):
        n = f.get('name')
        names.append(n)
        out[n] = {'type': f.get('type'),
                  'multi': f.get('multiValued') == 'true'}
    return out, names


def walk(doc, core, seen):
    """Collect {field: set(python types)} and whether it ever arrives as a list."""
    for k, v in doc.items():
        if isinstance(v, list) and v and isinstance(v[0], dict):
            if k in NEST_CONTAINERS[core]:
                for child in v:
                    walk(child, core, seen)
                continue
        rec = seen.setdefault(k, {'types': set(), 'listed': False, 'n': 0})
        rec['n'] += 1
        if isinstance(v, list):
            rec['listed'] = True
            for item in v:
                rec['types'].add(type(item))
        else:
            rec['types'].add(type(v))


def main():
    if '--compile' in sys.argv:
        print('Recompiling ...')
        subprocess.run([sys.executable, 'Compile_Biochemistry_for_SOLR.py'],
                       cwd=os.path.join(SOLR, 'compilation'), check=True)

    failures = []
    for core, docs_path in CORES.items():
        print(f'\n── {core} ──')
        if not os.path.exists(docs_path):
            failures.append(f'{core}: {docs_path} missing; run with --compile')
            print('  FAIL — compiler output missing')
            continue
        age = (time.time() - os.path.getmtime(docs_path)) / 3600.0
        print(f'  compiler output is {age:.1f} h old')

        # 1 + 2: schema parses, no duplicate declarations
        try:
            decl, names = declared_fields(core)
        except ET.ParseError as e:
            failures.append(f'{core}: schema.xml does not parse: {e}')
            print('  FAIL — schema.xml does not parse')
            continue
        dupes = sorted({n for n in names if names.count(n) > 1})
        if dupes:
            failures.append(f'{core}: duplicate field declarations {dupes}')
        print(f'  {len(decl)} fields declared, {len(dupes)} duplicated')

        docs = json.load(open(docs_path))
        seen = {}
        ids, dup_ids, childless = set(), [], 0
        for d in docs:
            walk(d, core, seen)
            stack = [(d, True)]
            while stack:
                cur, is_parent = stack.pop()
                i = cur.get('id')
                if not i:
                    childless += 1
                elif i in ids:
                    dup_ids.append(i)
                else:
                    ids.add(i)
                if not is_parent and not cur.get('doc_type'):
                    childless += 1
                for k, v in cur.items():
                    if k in NEST_CONTAINERS[core] and isinstance(v, list):
                        stack.extend((c, False) for c in v)

        # 3: undeclared
        undeclared = sorted(k for k in seen
                            if k not in decl and k not in NEST_CONTAINERS[core])
        if undeclared:
            for k in undeclared:
                failures.append(f'{core}: field {k!r} emitted on {seen[k]["n"]:,} '
                                f'docs but not declared in schema.xml')
        print(f'  {len(seen)} fields emitted, {len(undeclared)} undeclared')

        # whitelist discipline: a container must actually be used as a container
        for c in NEST_CONTAINERS[core]:
            if c in decl:
                failures.append(f'{core}: {c!r} is whitelisted as a nest container '
                                f'but is also declared as a field')

        # 4 + 5: types and cardinality
        bad_type = bad_card = 0
        for k, rec in seen.items():
            spec = decl.get(k)
            if not spec:
                continue
            want = NUMERIC.get(spec['type'])
            if want:
                for t in rec['types']:
                    if t is bool or not issubclass(t, (int, float)):
                        failures.append(f'{core}: {k!r} declared {spec["type"]} '
                                        f'but emitted {t.__name__}')
                        bad_type += 1
                    elif want is int and t is float:
                        failures.append(f'{core}: {k!r} declared int but emitted float')
                        bad_type += 1
            if spec['type'] == 'boolean':
                for t in rec['types']:
                    if t is not bool:
                        failures.append(f'{core}: {k!r} declared boolean but '
                                        f'emitted {t.__name__}')
                        bad_type += 1
            if rec['listed'] and not spec['multi']:
                failures.append(f'{core}: {k!r} emitted as a list but not '
                                f'declared multiValued')
                bad_card += 1
        print(f'  type mismatches {bad_type}, cardinality mismatches {bad_card}')

        # 6: key integrity
        if dup_ids:
            failures.append(f'{core}: {len(dup_ids)} duplicate ids, '
                            f'e.g. {dup_ids[:3]}')
        if childless:
            failures.append(f'{core}: {childless} docs missing id or doc_type')
        print(f'  {len(ids):,} unique ids, {len(dup_ids)} duplicates, '
              f'{childless} missing id/doc_type')

        # 7: refusal sentinels must not be indexed as energies
        def numeric(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                return None
        flat_sent = child_sent = 0
        for d in docs:
            v = numeric(d.get('deltag'))
            if v is not None and abs(v) >= SENTINEL_DG:
                flat_sent += 1
            for c in d.get('thermodynamics', []) or []:
                e, s = numeric(c.get('energy')), numeric(c.get('error'))
                if (e is not None and abs(e) >= SENTINEL_DG) or \
                   (s is not None and abs(s) >= SENTINEL_SIGMA):
                    child_sent += 1
        # Reported, not failed -- indexing sentinels is a deliberate choice.
        print(f'  refusal sentinels indexed: {flat_sent:,} parent, '
              f'{child_sent:,} child   (by decision, not a failure)')

    print()
    if failures:
        print(f'RESULT: FAIL — {len(failures)} issue(s)')
        for f in failures[:25]:
            print(f'  {f}')
        if len(failures) > 25:
            print(f'  ... and {len(failures) - 25} more')
        sys.exit(1)
    print('RESULT: PASS')
    sys.exit(0)


if __name__ == '__main__':
    main()
