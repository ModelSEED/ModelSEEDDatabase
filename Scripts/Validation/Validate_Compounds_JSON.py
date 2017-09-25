import Scripts.Validation.Schemas as schemas
from Scripts.Validation.error_reporting import find_new_errors, report_errors, get_master_errors
from repostat.stash import StatStash
from jsonschema import Draft4Validator
from collections import defaultdict
import argparse
import json


def validate_schema(_compounds, verbose):
    v = Draft4Validator(schemas.compounds)
    schema_errors = sorted(v.iter_errors(_compounds), key=lambda e: e.path)
    if verbose:
        for e in schema_errors:
            print("Schema Validation ERROR at compounds[{}]: {}"
                  .format("][".join(e.absolute_path), e.message))
    return len(schema_errors)


def check_dups(_compounds, verbose, unique_fields=('id', 'abbreviation',
                                                   'name', 'structure')):
    unique_values = dict([(x, defaultdict(list)) for x in unique_fields])
    for id, comp in _compounds.items():
        if comp['is_obsolete'] == "1":
            continue
        for key in unique_values:
            unique_values[key][comp[key]].append(id)

    duplicated = defaultdict(int)
    for value_type in unique_values:
        for key, ids in unique_values[value_type].items():
            if len(ids) > 1 and key != 'null':
                if verbose:
                    print("Duplicate {}: {} in {}".format(value_type, key, ids))
                duplicated["duplicated_"+value_type] += len(ids)-1
    return duplicated


if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser()
    parser.add_argument('cpdfile', help='path to compounds file', action='store')
    parser.add_argument('-d', '--print-duplicates', help='Print all duplicated values', action='store_true', dest='print_dups', default=False)
    parser.add_argument('-s', '--print-schema', help='Print all schema violations', action='store_true', dest='print_schema', default=False)
    args = parser.parse_args()
    compounds = json.load(open(args.cpdfile))
    errors = dict()
    errors['schema_violations'] = validate_schema(compounds, args.print_schema)
    errors.update(check_dups(compounds, args.print_dups))
    print("\nError counts")
    for error_type, count in errors.items():
        print("\t{}: {}".format(error_type, count))
    new_errors = find_new_errors('compounds', errors)
    report_errors('compounds', errors)
    print("New errors detected: " + ", ".join(new_errors))
    exit(len(new_errors))
