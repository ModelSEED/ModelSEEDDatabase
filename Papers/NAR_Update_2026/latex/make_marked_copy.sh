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
# the journal page limit). Per NAR's own instruction ("Type your revised
# text in red font, please DO NOT use red highlight") the marked copy shows
# ONLY additions, in plain red with no underline -- deletions are
# completely omitted (--no-del below), not struck-through, so there is no
# blue text anywhere in the output.
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

latexdiff --flatten --no-del \
  "$WORK/old/Papers/NAR_Update_2026/latex/$DOC.tex" \
  "$WORK/new/Papers/NAR_Update_2026/latex/$DOC.tex" \
  > "$HERE/${DOC}_diff.tex"

# NAR: "Type your revised text in red font, please DO NOT use red
# highlight" -- plain red text, no underline/wave/strikethrough decoration.
# --no-del above already drops deleted text entirely (no blue anywhere).
# Swap latexdiff's default \DIFadd (blue, wavy-underlined via \uwave) for
# plain \color{red} with no decoration.
python3 - "$HERE/${DOC}_diff.tex" <<'PYEOF'
import re, sys
path = sys.argv[1]
with open(path) as f:
    text = f.read()

# latexdiff emits \DIFadd directly in most document classes, but for classes
# it detects as hyperref-wrapped (e.g. this repo's supplementary.tex, which
# loads hyperref before latexdiff's preamble) it instead defines \DIFaddtex
# and routes \DIFadd through \texorpdfstring. Try the direct macro first,
# then the *tex variant, so this one script handles both main.tex and
# supplementary.tex.
pairs = [
    (r"\providecommand{\DIFadd}[1]{{\protect\color{blue}\uwave{#1}}} %DIF PREAMBLE",
     r"\providecommand{\DIFadd}[1]{{\protect\color{red}{#1}}} %DIF PREAMBLE"),
    (r"\providecommand{\DIFaddtex}[1]{{\protect\color{blue}\uwave{#1}}} %DIF PREAMBLE",
     r"\providecommand{\DIFaddtex}[1]{{\protect\color{red}{#1}}} %DIF PREAMBLE"),
]

applied = False
for add_old, add_new in pairs:
    if add_old in text:
        text = text.replace(add_old, add_new)
        applied = True
        break

if not applied:
    sys.exit("make_marked_copy.sh: expected DIFadd preamble line not found -- "
             "latexdiff output format may have changed, update the swap in this script.")

# Safety net: --no-del is supposed to drop every deleted block, but a
# deletion word-adjacent to an addition on the same line (e.g. a title-block
# edit mixing \DIFdelbegin...\DIFdelend with \DIFaddbegin...\DIFaddend) has
# been observed to survive it. Scoped to the body only (after \begin{document})
# so the preamble's own \providecommand{\DIFdel}... definitions are never
# touched. \DIFdel{...} is removed with a brace-depth counter, not a regex --
# its argument can itself contain braced groups (\textbf{...} etc.), which a
# naive regex mishandles and corrupts the file.
def strip_body_deletions(body):
    # \DIFdelbegin ... \DIFdelend spans, non-nesting (latexdiff never nests
    # these), so a simple non-greedy match across the whole body is safe.
    body = re.sub(r"\\DIFdelbegin\b.*?\\DIFdelend\b", "", body, flags=re.S)

    # Bare \DIFdel{...} calls (left outside a begin/end span): find each
    # occurrence and consume balanced braces by hand.
    out = []
    i = 0
    marker = r"\DIFdel{"
    while True:
        j = body.find(marker, i)
        if j == -1:
            out.append(body[i:])
            break
        out.append(body[i:j])
        k = j + len(marker)
        depth = 1
        while depth > 0:
            if body[k] == "{":
                depth += 1
            elif body[k] == "}":
                depth -= 1
            k += 1
        i = k  # drop body[j:k] entirely -- the whole \DIFdel{...} call
    return "".join(out)

# Match \begin{document} only at the start of a (whitespace-trimmed) line, so
# a comment mentioning "\begin{document}" earlier in the preamble (main.tex
# has one, documenting \draftmodefalse) is not mistaken for the real one --
# that mistake fed the preamble's own \DIFdel macro definitions into the
# stripper above and corrupted them.
m = re.search(r"^[ \t]*\\begin\{document\}", text, flags=re.M)
if not m:
    sys.exit("make_marked_copy.sh: \\begin{document} not found -- "
             "refusing to touch the preamble's \\DIFdel definitions.")
doc_start = m.start()
text = text[:doc_start] + strip_body_deletions(text[doc_start:])

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
