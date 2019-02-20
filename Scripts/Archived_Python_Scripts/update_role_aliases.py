import csv
import os
import sys


def make_rename_dict(rewrite_path):
    rename_mapping = {}
    for line in open(rewrite_path, encoding='utf-16'):
        try:
            k, v = line.split(" was replaced by ")
            rename_mapping[k[5:]] = v.strip()
        except ValueError:
            continue
    print('{} roles renamed'.format(len(rename_mapping)))
    return rename_mapping


def update_roles(role_path, rename_dict, new_path=None):
    if not new_path:
        new_path = role_path
    updated = 0
    tmp = 'roles.tmp'
    reader = csv.DictReader(open(role_path), dialect='excel-tab')
    writer = csv.DictWriter(open(tmp, 'w'), reader.fieldnames, dialect='excel-tab')
    writer.writeheader()
    for row in reader:
        if row['name'] in rename_dict:
            updated += 1
            if row['aliases'] == 'null':
                row['aliases'] = rename_dict[row['name']]
            else:
                row['aliases'] += ";" + rename_dict[row['name']]
        writer.writerow(row)
    print("{} role aliases updated".format(updated))
    os.rename(tmp, new_path)


if __name__ == "__main__":
    roles_path = os.path.dirname(sys.argv[0])+"/../../Annotations/Roles.tsv"
    role_rewrite_path = sys.argv[1]
    rename_dict = make_rename_dict(role_rewrite_path)
    update_roles(roles_path, rename_dict)
