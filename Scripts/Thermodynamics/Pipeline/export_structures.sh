#!/usr/bin/env bash
# Export ModelSEED structure + pKa data from THIS repository at a pinned commit
# into a local snapshot, so the eQuilibrator cache build is reproducible and
# independent of whatever branch the working tree happens to be on.
#
# Usage:  export_structures.sh [ref]     (default: HEAD of the current checkout)
set -euo pipefail

SRC_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
REF=${1:-HEAD}
OUT="${EQUILIBRATOR_DIR:-/scratch/seaver/Claude_Projects/eQuilibrator}/data/msd_structures"

cd "$SRC_REPO"
SHA=$(git rev-parse "$REF")
mkdir -p "$OUT"

for s in KEGG ChEBI MetaCyc Rhea; do
    mkdir -p "$OUT/$s"
    # pKa ladders: export every marvin_*.tsv present (load_marvin_pkas merges
    # them, newest version winning per id). Protonation bundles: the consumed
    # one per sources.yaml (marvin_26.1_ph7.tsv today).
    mkdir -p "$OUT/$s/pkas" "$OUT/$s/protonations"
    git ls-tree --name-only "$SHA" -- "Biochemistry/Structures/$s/pkas/" | while read -r f; do
        git show "$SHA:$f" > "$OUT/$s/pkas/$(basename "$f")"
    done
    git ls-tree --name-only "$SHA" -- "Biochemistry/Structures/$s/protonations/" | while read -r f; do
        git show "$SHA:$f" > "$OUT/$s/protonations/$(basename "$f")"
    done
    for f in smiles.tsv inchikey.tsv; do
        git show "$SHA:Biochemistry/Structures/$s/$f" > "$OUT/$s/$f"
    done
done

mkdir -p "$OUT/ModelSEED/pkas"
git ls-tree --name-only "$SHA" -- "Biochemistry/Structures/ModelSEED/pkas/" | while read -r f; do
    git show "$SHA:$f" > "$OUT/ModelSEED/pkas/$(basename "$f")"
done
git show "$SHA:Biochemistry/Structures/sources.yaml" > "$OUT/sources.yaml"
git show "$SHA:Biochemistry/Structures/Unique_ModelSEED_Structures.txt" \
    > "$OUT/Unique_ModelSEED_Structures.txt"

cat > "$OUT/PROVENANCE.txt" <<PROV
Exported from $SRC_REPO
ref:    $REF
commit: $SHA
$(git log -1 --format='date:   %ci%nsubject: %s' "$SHA")
PROV

echo "exported structures at $SHA ($REF) to $OUT"
find "$OUT" -name '*.tsv' -o -name '*.txt' -o -name '*.yaml' | sort | while read -r f; do
    printf '  %-58s %8d lines\n' "${f#$OUT/}" "$(wc -l < "$f")"
done
