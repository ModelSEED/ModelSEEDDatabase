"""Validates Annotations and Media files"""

from csv import DictReader
import os
import argparse


def check_id_set(tsv_path):
    with open(tsv_path) as infile:
        unique, duplicates = set(), []
        for line in infile:
            val = line.split('\t')[0]
            if val in unique:
                duplicates.append(val)
            else:
                unique.add(val)
        return unique, duplicates


def validate_media_compounds(path, comp_set):
    media_compounds = set([x.strip() for x in open(path)])
    return media_compounds - comp_set


if __name__ == '__main__':
    script_dir = os.path.dirname(__file__)
    comp_ids, dup_comps = check_id_set(script_dir+"/../../Biochemistry/compounds.tsv")
    complex_ids, dup_complex = check_id_set(script_dir+"/../../Annotations/Complexes.tsv")
    role_ids, dup_roles = check_id_set(script_dir+"/../../Annotations/Roles.tsv")
    undefined_comps = validate_media_compounds(script_dir+"/../../Media/KBaseMedia.cpd",
                                               comp_ids)

    if undefined_comps:
        print("ERROR-Undefined Compounds: " + ", ".join(undefined_comps))
    if dup_complex:
        print("ERROR-Duplicated Complexes: " + ", ".join(dup_complex))
    if dup_roles:
        print("ERROR-Duplicated Roles: " + ", ".join(dup_roles))
    exit(any([undefined_comps, dup_roles, dup_complex]))