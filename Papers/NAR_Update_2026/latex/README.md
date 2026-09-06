# NAR 2026 update — LaTeX source

LaTeX build of `../MANUSCRIPT.md`, targeting **Nucleic Acids Research**,
Database Issue ([author guidance](https://academic.oup.com/nar/pages/Ms_Prep_Database)).

> **Journal requirements — page budget, figure formats and dpi, graphical
> abstract, reference style — are documented with citations in
> [`NAR_REQUIREMENTS.md`](NAR_REQUIREMENTS.md). Read that before changing
> formatting or preparing figures.**

## Layout

```
latex/
├── main.tex                     # ONLY file with journal formatting; \input's everything
├── references.bib               # numeric (NAR) style; entries marked [VERIFY] need checking
├── oup-authoring-template.cls   # vendored from CTAN v1.5 (LPPL) so this builds anywhere
├── oup-plain.bst                # numeric bibliography style (NAR)
├── oup-abbrvnat.bst             # author-year alternative, unused
├── figures/                     # figure sources; empty until the direction-study figure lands
└── sections/                    # one file per manuscript section — all prose lives here
```

`main.tex` carries the document class, the NAR requirements checklist, the draft
markers, and the ordered list of `\input` lines. **Do not put prose in
`main.tex`**; add or reorder `\input` lines instead.

NAR specifies the OUP template's **Modern Large** design, set as
`\documentclass[unnumsec,webpdf,modern,large]{oup-authoring-template}`. That
gives 9bp body text on 11.5pt leading, a **sans-serif** body (modern is the only
one of the three OUP designs that forces sans), 7.5bp footnotes, 210×276 mm
paper, and two columns. Don't change `modern` or `large` — see
[`NAR_REQUIREMENTS.md` §7](NAR_REQUIREMENTS.md#required-design-modern-large).

### Sections

| File | Manuscript section |
|---|---|
| `abstract.tex` | Abstract (single paragraph — it is `\input` into `\abstract{}`) |
| `introduction.tex` | Introduction |
| `methods_collation.tex` … `methods_distribution.tex` | Materials and Methods (13 subsections, in reading order) |
| `results_growth.tex` … `results_connectivity_ontology.tex` | Results (7 subsections) |
| `discussion.tex` | Discussion |
| `data_availability.tex` | Data Availability |
| `author_contributions.tex` | Author contributions / Funding / Acknowledgements / COI |

Each file opens with a comment naming its source in `MANUSCRIPT.md` and the
originating file under `../drafts/`, so prose can be traced back.

## Building

```bash
latexmk -pdf main.tex
# or: pdflatex main && bibtex main && pdflatex main && pdflatex main
```

No TeX installation is required to edit — the project is a plain Overleaf
upload (the `.cls` is vendored, so nothing needs installing there either).

## Draft markers

The manuscript carries 78 `[TBD]` and 15 `[DRAFT PENDING]` placeholders. These
map to three macros defined in `main.tex`:

- `\TBD{}` / `\TBD[note]` — red, inline
- `\DRAFTPENDING{what is missing}` — orange, section not yet written
- `\NUMBERSPENDING{what is missing}` — orange, prose ready but numbers pending

Set `\draftmodefalse` in `main.tex` to make every marker disappear for a clean
submission PDF. Grep for remaining work with:

```bash
grep -rn 'TBD\|DRAFTPENDING\|NUMBERSPENDING' sections/
```

## Not carried over

The "Manuscript-shape health check" table at the end of `MANUSCRIPT.md` is a
project-tracking artifact, not manuscript content, so it is deliberately not in
the LaTeX. It stays in `MANUSCRIPT.md`.

## Before submission

Full sourced checklist: [`NAR_REQUIREMENTS.md`](NAR_REQUIREMENTS.md) (also
mirrored as a short list at the top of `main.tex`). The items needing action
beyond filling `[TBD]`s:

1. **Title** — NAR wants the database name as the first word; the title
   currently opens with "The".
2. **Graphical abstract** — mandatory, 5:2 ratio, min 127×50 mm, 300–600 dpi,
   original artwork. Not embedded here; submitted as a separate file.
3. **Length** — update papers run 4–6 typeset pages. This draft is over by
   design; trim once the `[TBD]`s resolve.
4. **References** — resolve every `[VERIFY]` note in `references.bib`.
5. **Six suggested referees** — supplied in ScholarOne, not in this source.
