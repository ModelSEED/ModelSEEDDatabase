#!/usr/bin/env bash
# Regenerate the referee-facing "marked" (tracked-changes) copy of the
# manuscript using latexdiff, word-level, against a chosen base ref.
#
# Usage (from this directory, i.e. Papers/NAR_Update_2026/latex/):
#   ./make_marked_copy.sh [BASE_REF] [DOC]
# BASE_REF defaults to origin/dev. DOC selects which top-level document to
# diff and defaults to "main" (main.tex -> main_diff.tex/pdf). Pass
# "supplementary" to build the supplement's redline instead
# (supplementary.tex -> supplementary_diff.tex/pdf) -- needed whenever a PR
# touches supplement/*.tex, since main.tex and supplementary.tex are two
# separate compiled documents that each \input their own sections.
#
# What it does:
#   1. Exports the full Papers/NAR_Update_2026 tree at BASE_REF ("old") and
#      takes the current working tree as "new".
#   2. Runs `latexdiff --flatten old/DOC.tex new/DOC.tex`, which expands
#      every \input (sections/*.tex, supplement/*.tex) into one file and
#      injects the DIFadd/DIFdel preamble macros -- this is why we diff the
#      whole assembled document rather than section files individually.
#   3. Drops the result as DOC_diff.tex next to DOC.tex (so its
#      ../figures/... paths still resolve) and compiles it:
#      pdflatex -> bibtex -> pdflatex -> pdflatex.
#
# Requires latexdiff (TeX Live 2026: /scratch/seaver/texlive/2026/bin/x86_64-linux).
# Output: DOC_diff.pdf (word-level track-changes copy, NOT counted against
# the journal page limit -- deletions still take up visual space as
# strikethrough text, so this PDF is typically longer than either version).
#
# The plain submission copy is unaffected by this script: build it with
# `pdflatex main_clean.tex` (see README.md) as before.

set -euo pipefail

BASE_REF="${1:-origin/dev}"
DOC="${2:-main}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../../.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "Base ref: $BASE_REF"
echo "Document: $DOC.tex"
echo "Scratch:  $WORK"

mkdir -p "$WORK/old" "$WORK/new/Papers/NAR_Update_2026"
( cd "$REPO_ROOT" && git archive "$BASE_REF" -- Papers/NAR_Update_2026 ) | tar -x -C "$WORK/old"
cp -r "$REPO_ROOT/Papers/NAR_Update_2026/latex" "$WORK/new/Papers/NAR_Update_2026/"
cp -r "$REPO_ROOT/Papers/NAR_Update_2026/figures" "$WORK/new/Papers/NAR_Update_2026/"
find "$WORK/new/Papers/NAR_Update_2026/latex" \
  \( -name '*.aux' -o -name '*.bbl' -o -name '*.blg' -o -name '*.log' -o -name '*.out' \) -delete

export PATH="/scratch/seaver/texlive/2026/bin/x86_64-linux:$PATH"

latexdiff --flatten \
  "$WORK/old/Papers/NAR_Update_2026/latex/$DOC.tex" \
  "$WORK/new/Papers/NAR_Update_2026/latex/$DOC.tex" \
  > "$HERE/${DOC}_diff.tex"

# NAR requires that text changed in response to referee comments be shown in
# red. latexdiff's default UNDERLINE style colors additions blue and
# deletions red -- the opposite of what NAR wants. Swap the two auto-generated
# preamble macros in place: additions (\DIFadd) become red, deletions
# (\DIFdel) stay struck-through but become blue.
python3 - "$HERE/${DOC}_diff.tex" <<'PYEOF'
import re, sys
path = sys.argv[1]
with open(path) as f:
    text = f.read()

# latexdiff emits \DIFadd/\DIFdel directly in most document classes, but for
# classes it detects as hyperref-wrapped (e.g. this repo's supplementary.tex,
# which loads hyperref before latexdiff's preamble) it instead defines
# \DIFaddtex/\DIFdeltex and routes \DIFadd/\DIFdel through \texorpdfstring.
# Try the direct pair first, then the *tex pair, so this one script handles
# both main.tex and supplementary.tex.
pairs = [
    (r"\providecommand{\DIFadd}[1]{{\protect\color{blue}\uwave{#1}}} %DIF PREAMBLE",
     r"\providecommand{\DIFadd}[1]{{\protect\color{red}\uwave{#1}}} %DIF PREAMBLE",
     r"\providecommand{\DIFdel}[1]{{\protect\color{red}\sout{#1}}} %DIF PREAMBLE",
     r"\providecommand{\DIFdel}[1]{{\protect\color{blue}\sout{#1}}} %DIF PREAMBLE"),
    (r"\providecommand{\DIFaddtex}[1]{{\protect\color{blue}\uwave{#1}}} %DIF PREAMBLE",
     r"\providecommand{\DIFaddtex}[1]{{\protect\color{red}\uwave{#1}}} %DIF PREAMBLE",
     r"\providecommand{\DIFdeltex}[1]{{\protect\color{red}\sout{#1}}} %DIF PREAMBLE",
     r"\providecommand{\DIFdeltex}[1]{{\protect\color{blue}\sout{#1}}} %DIF PREAMBLE"),
]

applied = False
for add_old, add_new, del_old, del_new in pairs:
    if add_old in text and del_old in text:
        text = text.replace(add_old, add_new).replace(del_old, del_new)
        applied = True
        break

if not applied:
    sys.exit("make_marked_copy.sh: expected DIFadd/DIFdel preamble lines not found -- "
             "latexdiff output format may have changed, update the swap in this script.")

with open(path, "w") as f:
    f.write(text)
PYEOF

cd "$HERE"
pdflatex -interaction=nonstopmode "${DOC}_diff.tex" >/dev/null
bibtex "${DOC}_diff" >/dev/null || true
pdflatex -interaction=nonstopmode "${DOC}_diff.tex" >/dev/null
pdflatex -interaction=nonstopmode "${DOC}_diff.tex"

echo
echo "Wrote $HERE/${DOC}_diff.tex and ${DOC}_diff.pdf"
grep -o "Output written on ${DOC}_diff.pdf.*pages.*" "${DOC}_diff.log" || true
