"""Generate before/after structure images for applied corrections."""
import csv
import os
import sys
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CORRECTIONS_LOG = os.path.join(
    os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', '..')),
    'Cleanup', 'pubchem_corrections_log.tsv')
IMAGES_DIR = os.path.join(SCRIPT_DIR, 'struct_imgs')

sys.path.insert(0, os.path.join(SCRIPT_DIR, '..', '..', 'Libs', 'Python'))
from BiochemPy import Compounds

COMPOUND_FILES = sorted(
    f for f in os.listdir(os.path.join(SCRIPT_DIR, '..', '..', 'Biochemistry'))
    if f.startswith('compound_') and f.endswith('.tsv'))
DB_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', 'Biochemistry'))


def load_names():
    names = {}
    for fn in COMPOUND_FILES:
        with open(os.path.join(DB_ROOT, fn)) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                names[row['id']] = row.get('name', '')
    return names


def load_corrections():
    corrections = defaultdict(dict)
    with open(CORRECTIONS_LOG) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            cpd_id = row['cpd_id']
            field = row['field']
            if field == 'SMILE':
                corrections[cpd_id]['old_smiles'] = row['old_value']
                corrections[cpd_id]['new_smiles'] = row['new_value']
                corrections[cpd_id]['result_type'] = row['result_type']
                corrections[cpd_id]['strategy'] = row['strategy']
                corrections[cpd_id]['pubchem_cid'] = row['pubchem_cid']
            elif field == 'Formula/Charge':
                corrections[cpd_id]['old_fc'] = row['old_value']
                corrections[cpd_id]['new_fc'] = row['new_value']
                if 'result_type' not in corrections[cpd_id]:
                    corrections[cpd_id]['result_type'] = row['result_type']
    return dict(corrections)


def render_correction(cpd_id, corr, name, font_large, font_small):
    old_smi = corr.get('old_smiles', '')
    new_smi = corr.get('new_smiles', '')
    result_type = corr.get('result_type', '')

    mol_old = Chem.MolFromSmiles(old_smi) if old_smi else None
    mol_new = Chem.MolFromSmiles(new_smi) if new_smi else None
    if mol_old is None and mol_new is None:
        return False

    mol_size = (350, 300)
    img_old = (Draw.MolToImage(mol_old, size=mol_size) if mol_old
               else Image.new("RGB", mol_size, "white"))
    img_new = (Draw.MolToImage(mol_new, size=mol_size) if mol_new
               else Image.new("RGB", mol_size, "white"))

    margin, header_h, gap = 10, 50, 20
    w = mol_size[0] * 2 + gap + margin * 2
    footer_h = 55
    h = mol_size[1] + header_h + footer_h + margin * 2
    canvas = Image.new("RGB", (w, h), "white")
    draw = ImageDraw.Draw(canvas)

    label = result_type.replace('_', ' ').title()
    draw.text((margin, 5), f"{cpd_id}  [{result_type}]",
              fill="black", font=font_large)
    draw.text((w - margin - 80, 5), "CORRECTED",
              fill="green", font=font_large)

    left_x = margin
    right_x = margin + mol_size[0] + gap
    draw.text((left_x + mol_size[0] // 2 - 25, header_h - 18),
              "Before", fill="blue", font=font_large)
    draw.text((right_x + mol_size[0] // 2 - 20, header_h - 18),
              "After", fill=(34, 139, 34), font=font_large)

    canvas.paste(img_old, (left_x, header_h))
    canvas.paste(img_new, (right_x, header_h))

    draw.line([(left_x + mol_size[0] + gap // 2, header_h),
               (left_x + mol_size[0] + gap // 2, header_h + mol_size[1])],
              fill="gray", width=1)

    footer_y = header_h + mol_size[1] + 5
    info = f"Name: {name[:60]}"
    draw.text((margin, footer_y), info, fill="black", font=font_small)
    strategy = corr.get('strategy', '')
    cid = corr.get('pubchem_cid', '')
    info2 = f"Strategy: {strategy}   CID: {cid}"
    draw.text((margin, footer_y + 16), info2, fill="black", font=font_small)
    if 'old_fc' in corr:
        info3 = f"Formula/Charge: {corr['old_fc']} -> {corr['new_fc']}"
        draw.text((margin, footer_y + 32), info3,
                  fill="black", font=font_small)

    prefix = {'STEREO_DIFF': 'stereo', 'PKA_CORRECTION': 'pka',
              'PROTONATION_DIFF_CORRECTED': 'prot',
              'CONSISTENCY_FIX': 'consistency'}.get(result_type, 'other')
    fname = f"{prefix}_corrected_{cpd_id}.png"
    canvas.save(os.path.join(IMAGES_DIR, fname))
    return True


def main():
    max_per_type = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    corrections = load_corrections()
    names = load_names()
    print(f"Loaded {len(corrections)} corrections")

    try:
        font_large = ImageFont.truetype(
            "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", 16)
        font_small = ImageFont.truetype(
            "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", 12)
    except (OSError, IOError):
        font_large = font_small = ImageFont.load_default()

    by_type = defaultdict(list)
    for cpd_id, corr in corrections.items():
        if 'old_smiles' in corr:
            by_type[corr['result_type']].append(cpd_id)

    os.makedirs(IMAGES_DIR, exist_ok=True)

    for result_type, cpd_ids in sorted(by_type.items()):
        count = 0
        for cpd_id in sorted(cpd_ids):
            if count >= max_per_type:
                break
            name = names.get(cpd_id, '')
            if render_correction(cpd_id, corrections[cpd_id], name,
                                 font_large, font_small):
                count += 1
                print(f"  {result_type}: {cpd_id} ({name[:40]})")
        print(f"{result_type}: {count} images")


if __name__ == "__main__":
    main()
