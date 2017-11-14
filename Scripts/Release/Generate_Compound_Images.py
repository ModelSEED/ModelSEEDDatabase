import os
import csv
import platform
import argparse
from subprocess import call
from shutil import move


mol_convert_path = "/Applications/MarvinSuite/bin/molconvert"

def generate_image_files(compounds_file, path, structure_key='structure', dir_depth=0,
                         img_type='svg', formating='-a,nosource,w500,h500'):
    """Generates image files for compounds in database using ChemAxon's
    MolConvert.

    :param compounds_file: File containing compounds to be converted
    :type compounds_file: str
    :param path: Target directory for image files
    :type path: str
    :param query: A query to limit the number of files generated
    :type query: dict
    :param dir_depth: The number of directory levels to split the
        compounds into for files system efficiency. Ranges from 0 (all in
        top level directory to the length of the compounds_file name (40 for MINE
        hashes)
    :type dir_depth: int
    :param img_type: The type of image compounds_file to be generated. See molconvert
        documentation for valid options
    :type img_type: str

    """
    ids = []
    structures = []
    tmp_structures = path[:-4] + 'structures.tmp'

    if not os.path.exists(path):
        os.mkdir(path)

    # build a new structure file with no missing structures
    with open(compounds_file) as infile:
        r = csv.DictReader(infile, dialect='excel-tab')
        for line in r:
            if not line[structure_key] or line[structure_key] == 'null':
                continue
            structures.append(line[structure_key])
            ids.append(line['id'])
    with open(tmp_structures, 'w') as outfile:
            outfile.write("\n".join(structures))
    print("Images for %s structures will be generated" % len(structures))
    if platform.system() == 'Windows':
        rc = call(['cmd', "/C molconvert -mo %s/.%s %s:%s %s"
                   % (path, img_type, img_type, formating, tmp_structures)],
                  shell=True)
    else:
        rc = call([mol_convert_path + " -mgo %s/.%s %s:%s %s" % (
            path, img_type, img_type, formating, tmp_structures)], shell=True)
    if rc:
        raise RuntimeError("molconvert returned %s" % rc)
    os.remove(tmp_structures)

    for i, _id in enumerate(ids):
        old = os.path.join(path, "%s.%s" % ((i + 1), img_type))
        new = path
        for j in range(0, dir_depth):
            new = os.path.join(new, _id[j])
        if not os.path.exists(new):
            os.makedirs(new)
        new = os.path.join(new, _id + '.' + img_type)
        if os.path.isfile(old):
            move(old, new)



def safe_generate_image_files(compounds_file, path, structure_key='structure',
                             dir_depth=0, img_type='svg',
                             formating='-a,nosource,w500,h500'):
    """This is a slower way of generateing image files because it uses a 
    subprocess call for each compound but this does make it more robust to 
    errors in the compound structures

    :param compounds_file: File containing compounds to be converted
    :type compounds_file: str
    :param path: Target directory for image files
    :type path: str
    :param query: A query to limit the number of files generated
    :type query: dict
    :param dir_depth: The number of directory levels to split the
        compounds into for files system efficiency. Ranges from 0 (all in
        top level directory to the length of the compounds_file name (40 for MINE
        hashes)
    :type dir_depth: int
    :param img_type: The type of image compounds_file to be generated. See molconvert
        documentation for valid options
    :type img_type: str

    """
    tmp_structures = os.path.dirname(__file__) + '/structures.tmp'
    path = os.path.dirname(__file__) + "/" + path
    if not os.path.exists(path):
        os.mkdir(path)
    with open(compounds_file) as infile:
        r = csv.DictReader(infile, dialect='excel-tab')
        for i,line in enumerate(r):
            if not line[structure_key] or line[structure_key] == 'null' or i < 600:
                continue
            with open(tmp_structures, 'w') as outfile:
                outfile.write(line[structure_key])
            _id = line['id'].replace('cpd','')
            new = path
            for j in range(0, dir_depth):
                new = os.path.join(new, _id[j])
            if not os.path.exists(new):
                os.makedirs(new)
            new = os.path.join(new, line['id'] + '.' + img_type)
            rc = call([mol_convert_path, '-Yo', new, img_type+":"+formating, tmp_structures])
    os.remove(tmp_structures)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate images for '
                                                 'compounds in a tsv')
    parser.add_argument('-i', '--inpath', type=str, default="../../Biochemistry/compounds.tsv",
                        help='The path to place the resulting output directory')
    parser.add_argument('-o', '--outpath', type=str, default="../../Objects/img",
                        help='The path to place the resulting output directory')
    parser.add_argument('-d', '--depth', type=int, default=0,
                        help='The desired depth of the directory tree')
    parser.add_argument('-t', '--type', type=str, default="png",
                        help='The desired output file type: png, svg or jpeg')
    parser.add_argument('-f', '--formatting', type=str, default="w500h500#00ffffff",
                        help='The dimensions of the output file')
    parser.add_argument('-k', '--kekulize', type=bool, default=True,
                        help='Should aromatic structures be represented in '
                             'Kekule form')
    args = parser.parse_args()
    format_string = args.formatting
    if args.kekulize:
        format_string += ',-a'
    safe_generate_image_files(args.inpath, args.outpath, dir_depth=args.depth,
                              img_type=args.type, formating=args.formatting)

