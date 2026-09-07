#!/bin/bash
# Regenerate the eQuilibrator energy tables from ModelSEED structures.
#
# This is the UPSTREAM stage. ../Rerun_Thermodynamics.sh ingests the tables this
# produces; it does not produce them. Run this only when the structures, the pKa
# layer or the training data have changed -- it needs the eQuilibrator working
# tree and takes roughly an hour.
#
#   EQUILIBRATOR_DIR=/path/to/eQuilibrator ./regenerate.sh
#
# Halts on the first failure. Each stage's output is the next stage's input, so
# continuing past an error would build on a partial artefact.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
: "${EQUILIBRATOR_DIR:?set EQUILIBRATOR_DIR to the eQuilibrator working tree}"
export EQUILIBRATOR_DIR
PY="${PYTHON:-python3}"

# -2..16, not 0..14: eQuilibrator counts sites above the reported pH to place
# the major microspecies, so a group at pKa 14.9 is protonated at pH 7. Without
# these, equilibrator-assets clips user-supplied pKas regardless of the caller.
export EQUILIBRATOR_PKA_MIN="${EQUILIBRATOR_PKA_MIN:--2}"
export EQUILIBRATOR_PKA_MAX="${EQUILIBRATOR_PKA_MAX:-16}"

# Pinned to the conditions the released tables report. pMg 3.0 rather than the
# standard preset's 10: free intracellular Mg2+ is of order 1 mM, and ATP is
# predominantly Mg-bound in vivo.
COND=(--p-h 7.0 --ionic-strength 0.25 --p-mg 3.0 --temperature 298.15)
CACHE="$EQUILIBRATOR_DIR/data/cache_final/compounds.sqlite"
PARAMS="$EQUILIBRATOR_DIR/data/cc_params_final.npz"

cd "$HERE"
echo "=== 1/8 seed mapping";        $PY build_seed_mapping.py
echo "=== 2/8 must-carry set";      $PY must_carry_over.py
echo "=== 3/8 resolve pKas";        $PY write_resolved_pkas.py
echo "=== 4/8 build cache";         $PY build_modelseed_cache.py --no-carry-over \
      --pkas "$EQUILIBRATOR_DIR/data/resolved_pkas.tsv" --out "$CACHE"
echo "=== 5/8 refit";               $PY train_path_b.py --cache "$CACHE" --out "$PARAMS"
echo "=== 6/8 reaction energies";   $PY generate_modelseed_energies.py \
      --cache "$CACHE" --params "$PARAMS" "${COND[@]}" \
      --out "$EQUILIBRATOR_DIR/data/modelseed_energies.tsv"
echo "=== 7/8 compound energies";   $PY generate_modelseed_compound_energies.py \
      --cache "$CACHE" --params "$PARAMS" "${COND[@]}" \
      --out "$EQUILIBRATOR_DIR/data/modelseed_formation_energies.tsv"

# Publish the reaction table, then ingest. update_modelseed_thermodynamics.py
# verifies both inputs were generated from $CACHE and refuses otherwise.
cp "$EQUILIBRATOR_DIR/data/modelseed_energies.tsv" \
   "$REPO/Biochemistry/Thermodynamics/eQuilibrator/ModelSEED_Reaction_Energies.tsv"
echo "=== 8/8 ingest into JSON";    $PY update_modelseed_thermodynamics.py --cache "$CACHE"

echo
echo "done. Verify with:  $PY verify_regen.py"
echo "Reactions newly covered and needing a reversibility call are listed in"
echo "  \$EQUILIBRATOR_DIR/data/newly_covered_needs_reversibility.tsv"
