# Working in this directory

LaTeX source for the ModelSEED 2026 update manuscript, targeting **Nucleic Acids
Research, Database Issue**.

## Hard rules — NAR will reject or mangle these

- **Never add a `\footnote`.** NAR explicitly prohibits footnotes. Use
  parenthetical text or a table note instead.
- **Never enable line numbering** (`lineno`, `\linenumbers`). Also prohibited.
- **Do not change the document class options.** NAR requires the OUP template's
  Modern Large design: `[unnumsec,webpdf,modern,large]`.
- **Do not add LaTeX packages** unless genuinely necessary — OUP asks authors to
  "refrain from adding further packages or macros if possible". The class
  already provides amsmath, graphicx, hyperref, natbib, xcolor, tikz, and its
  own `\toprule`/`\midrule`/`\botrule` (so never load booktabs — it clashes).
- **Do not put prose in `main.tex`.** Every section lives in `sections/`.

## Before changing anything

Read [`NAR_REQUIREMENTS.md`](NAR_REQUIREMENTS.md). It carries NAR's and OUP's
actual rules with citations and verbatim quotes — page budget, figure formats
and dpi, the mandatory graphical abstract, reference style, initial-submission
rules, and the submission checklist. It is the authority.

Two sections deserve particular attention:

- **§11** documents an unresolved conflict between OUP's `oup-plain.bst`
  (which sorts alphabetically) and NAR's order-of-appearance rule.
- **§12** lists ten formatting assumptions that are *guesses*, not sourced
  requirements, because NAR's Manuscript Preparation page is unreachable.
  Don't cite them back as if they were settled.

Do not restate requirements from memory, and do not add one to that file
without a source. If you learn a requirement, add it there with its citation.

Structure and build instructions: [`README.md`](README.md).
