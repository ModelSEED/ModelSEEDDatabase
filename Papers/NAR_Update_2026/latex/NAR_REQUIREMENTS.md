# NAR submission requirements — reference

Everything below is sourced. Each claim carries a link to where it came from and,
where the wording matters, a verbatim quote. **Do not add a requirement to this
file without a citation.**

Target: **Nucleic Acids Research, Database Issue — update paper.**
All pages checked **2026-08-26**.

---

## Quick answers

| Question | Answer | Source |
|---|---|---|
| How much text? | **4–6 typeset journal pages.** No word limit exists. | [DB Issue](https://academic.oup.com/nar/pages/Ms_Prep_Database) |
| How many figures? | **No stated limit.** The page budget is the real constraint. | [DB Issue](https://academic.oup.com/nar/pages/Ms_Prep_Database) |
| Figure format? | Raster → uncompressed **.tif**; vector → **.eps/.svg/.pdf** with embedded fonts. | [OUP artwork PDF](#sources) |
| Figure resolution? | **300 dpi** colour half-tone · **600** greyscale · **600–900** combination/line art · **1200** mono line art. | [OUP artwork PDF](#sources) |
| Graphical abstract? | **Mandatory.** 5:2, ≥127×50 mm, TIF/EPS/PDF, 300–600 dpi. | [DB Issue](https://academic.oup.com/nar/pages/Ms_Prep_Database) |
| LaTeX design? | OUP template, **'Modern Large'** — `[unnumsec,webpdf,modern,large]`. | [§7](#required-design-modern-large) |
| Initial submission? | **One .pdf**, figures/tables embedded. Supplementary separate. **No line numbers, no footnotes.** | [§10](#10-initial-submission) |
| Known problem? | **`oup-plain.bst` sorts alphabetically**, conflicting with NAR's order-of-appearance rule. | [§11](#11-known-conflict-bibliography-ordering) |
| Deadline? | **15 September** for update papers. | [DB Issue](https://academic.oup.com/nar/pages/Ms_Prep_Database) |

---

## 1. Length

Update papers:

> "should typically be no more than 4-6 journal pages in length"

New database papers (not our case):

> "New submissions to *NAR* are typically 4-5 typeset journal pages in length, but authors are urged to be succinct"

Source: [NAR Database Issue Guidelines](https://academic.oup.com/nar/pages/Ms_Prep_Database).
Contact the Executive Editor before submitting anything longer.

**No word count is specified on any NAR page.** NAR budgets in *typeset pages*.
The only way to know where we stand is to compile and count. The sibling
[Web Server Issue page](https://academic.oup.com/nar/pages/Submission_Webserver)
uses the same unit — "Submissions should typically be 4-5 printed journal pages
in length" — confirming this is how NAR expresses limits, not an oversight on
one page.

---

## 2. Figures — how many

**No numeric limit is stated** on the Database Issue page, the Web Server Issue
page, the [Author Guidelines](https://academic.oup.com/nar/pages/author-guidelines),
or the [Methods Guidelines](https://academic.oup.com/nar/pages/methods-guidelines).

Content *is* constrained. From the Database Issue page, the database home page

> "should not be used as a figure in the main text article to be typeset, but a representative screen dump"

of a query output is permitted.

---

## 3. Figures — format, resolution, colour

From OUP's *Guidance for preparing artwork* (see [Sources](#sources)).

**Read this caveat first**, it is the PDF's own opening line:

> "Below are tips rather than strict rules, as image content, source material and available software may constrain what you are able to achieve. Please use these guidelines alongside any specific instructions provide on the website of the journal to which you are submitting."

### File formats

- "Save raster images (photographs, scans) as uncompressed .tif to avoid quality loss"
- ".jpg/.png are acceptable for raster images but may be lower resolution than .tif"
- "Save vector images (diagrams, shapes, text) as .eps/.svg/.pdf and embed fonts"
- "For .svg, convert text to shapes/paths to ensure consistent browser display"
- "Avoid .bmp, .gif, and native application formats"
- "Save images created in MS Office with 'Print to PDF' to preserve quality and format"

### Resolution

Set the intended print size in your software so the resulting dpi is:

| Image type | Minimum dpi |
|---|---|
| Colour half-tones | "at least 300dpi" |
| Greyscale half-tones | "at least 600dpi" |
| Combination half-tones and line art | "600–900dpi" |
| Monochrome line art | "at least 1200dpi" |
| Pure vector | "no inherent resolution" |

- "Do not use 'up-sampling' in your image-editing software to artificially increase dpi"
- "Check your chosen journal's print size by using your PDF reader's measuring tool, to assess the maximum width of single and double-column images"
- "To reduce file size, crop white borders (minimum 2px) and flatten all layers"

### Colour

- "Either RGB or CMYK colour-space is acceptable"
- "RGB is recommended as it makes the best use of screen capabilities"
- "RGB images will be converted to CMYK for print"

### Fonts and lines

- "Use Arial, Times New Roman, Courier or Symbol fonts for accurate reproduction"
- "Ensure text is no less than 7pt"
- "Set line thickness between 0.25pt and 1pt"
- "Avoid pale colours such as yellow and colour combinations that may be difficult for colour-blind readers to distinguish (eg red–green), and favour bold contrasts"
- "Avoid using colour in isolation to convey meaning—consider textures, labels or additional text"

### File handling

- "Name figure files simply to match citation, eg fig1.tif, fig2.eps"
- "Provide multi-panel images in a single file"
- "Provide captions in the manuscript, not in the image file"
- "Cite all figures in sequence"
- "If combining multiple panels, avoid MS tools such as Powerpoint and use dedicated software such as Photoshop, GIMP, Illustrator or InkScape"

### At proof stage

- "PDF proofs reduce raster image resolution to 200dpi to manage online file size" — this is expected, not a defect in your figure
- "Vector images are preserved in the PDF proof to ensure no loss of quality"
- "Print size may be adjusted during typesetting to optimize layout"

**Practical consequence for this manuscript:** our figures are generated by
Python scripts, so export vector PDF/EPS rather than PNG wherever possible —
vector has no resolution ceiling, survives the proof stage intact, and sidesteps
the entire dpi table above. Matplotlib/Plotly `savefig(..., format='pdf')` with
fonts embedded. Use Arial to match the OUP list.

---

## 4. Graphical abstract — MANDATORY

> "Authors MUST provide a Graphical Abstract."

Specification, verbatim from the [Database Issue Guidelines](https://academic.oup.com/nar/pages/Ms_Prep_Database):

- "Size: 5:2 aspect ratio, 127x50mm or 5x2in minimum"
- "File Type: TIF, EPS or editable PDF"
- "Resolution: 300-600dpi minimum"
- "Orientation: landscape"
- "Font: Use a sans serif font such as Arial, 12–16 points"

Content requirements — it should "be simple", "be original i.e. not an existing
main or supplementary figure", "use colour", "use text sparingly, mainly for
labels", "read from top down or left to right", and must "not include trademarked
or copyrighted images or logos" (the example given: the text *UniProt* is fine,
"but not the logo").

It is submitted as a **separate file**, not embedded in the LaTeX.

---

## 5. Title, abstract, URL

- The database name "should ideally be the first word of the title."
  **Our title currently opens with "The"** — see the checklist in `main.tex`.
- A working URL must appear **in the abstract and in the article body**.
- Abstracts must stand alone: no citations to the reference list, no equations.

Source: [Database Issue Guidelines](https://academic.oup.com/nar/pages/Ms_Prep_Database).

---

## 6. References

- Cited in text by **sequential number in order of appearance**, listed
  numerically, in correct journal format.
- Excluded: items "submitted" or "in preparation", unpublished results, personal
  communications.
- Other databases: cite via their most recent published description. If none
  exists, "the URL goes in the body text rather than the reference list."
- **No maximum reference count is stated.**

Source: [Database Issue Guidelines](https://academic.oup.com/nar/pages/Ms_Prep_Database).

In this project the OUP class default (no `namedate`/`numbered` option) already
produces NAR-style numeric citations; the bibliography style is `oup-plain`.

---

## 7. Manuscript file formats and template

- Text "including references, figure legends and simple tables" may be `.pdf`,
  `.doc`, `.rtf`, or LaTeX. A PDF is acceptable at initial submission.
- NAR offers Word and LaTeX templates; the LaTeX one is the general OUP template.
- Submission is through ScholarOne Manuscripts.

Source: [Database Issue Guidelines](https://academic.oup.com/nar/pages/Ms_Prep_Database).

### Required design: Modern Large

> NAR recommends the OUP LaTeX template, available on Overleaf and as a
> downloadable package via OUP's *Preparing and submitting your manuscript*
> page. **Use the 'Modern Large' design.**

Source 7 (see [Sources](#sources)) — relayed by the corresponding author,
2026-08-26. Not independently retrieved; OUP's *Preparing and submitting your
manuscript* page returned navigation-only content on every automated fetch.

This is set in `main.tex` as:

```latex
\documentclass[unnumsec,webpdf,modern,large]{oup-authoring-template}
```

**Do not change `modern` or `large` without re-checking this section.**

What Modern Large produces, read from `oup-authoring-template.cls` v1.5:

| Property | Value | Class line |
|---|---|---|
| Body text | 9bp on 11.5pt leading | 175 |
| Body family | **sans-serif** — `modern` is the only design that forces `\sffamilyfont` | 174 |
| Footnotes | 7.5bp on 8 | 286 |
| Paper | 210 × 276 mm | 81–82 |
| Columns | two (`\twocolumn`) | 2521 |

For contrast, the template's own default (`contemporary, large`) is 8bp/11.5bp
with a serif body and 6.5bp footnotes. **Modern Large sets larger body type, so
it fits fewer words per page** — relevant against the 4–6 page budget in §1.

Main figures are single-column at this width; use `figure*` (and `table*`) to
span both columns.

We vendor [`oup-authoring-template`](https://ctan.org/pkg/oup-authoring-template)
v1.5 (2026-07-14, LPPL) so the project builds without installing anything.

The same class is published by OUP on Overleaf as
**[OUP General Template](https://www.overleaf.com/latex/templates/oup-general-template/ybpypwncdxyb)**
(source 11) — "A general template supporting journals published by Oxford
University Press." CTAN and Overleaf are the same template; we took the CTAN
copy because it can be vendored and version-pinned.

**Note when uploading to Overleaf:** because `oup-authoring-template.cls` is
committed here, an Overleaf project built from this directory uses *our* pinned
copy, not Overleaf's. That is deliberate — the build stays reproducible. If you
would rather track Overleaf's copy, delete our `.cls` and start from the gallery
template instead. The Overleaf listing showed a more recent update date
(early August 2026) than CTAN's v1.5 release (2026-07-14), so a small version
drift is possible; CTAN reported 1.5 as current when checked on 2026-08-26.

### Which template — and why not the legacy `NAR.cls`

A second, older template exists: a NAR-specific class used as
`\documentclass[a4,center,fleqn]{NAR}`. **We are not using it.** Its markup is
wholly incompatible with the OUP template:

| | legacy `NAR.cls` | `oup-authoring-template` (ours) |
|---|---|---|
| Class options | `a4, center, fleqn` | `unnumsec, webpdf, modern, large` |
| Volume / issue | `\jvolume`, `\jissue` | `\vol`, `\issue` |
| Authors | one `\author{}` block, manual `$^{1,*}$` marks | per-author `\author[n]{}` |
| Affiliations | one `\address{}` block | `\address[n]{}` per affiliation |
| Correspondence | `\footnote` inside `\author` | `\corresp[$\ast$]{}` |
| Dates | `\history{Received…; Revised…; Accepted…}` | `\received{}{}{}` etc. |
| Abstract | `\begin{abstract}…\end{abstract}` | `\abstract{…}` |
| Tables | `\tableparts{caption}{tabular}{note}`, `\colrule` | `\caption{}` + `\midrule` |
| Sections | numbered, ALL-CAPS headings | `unnumsec` |
| Bibliography | hand-written `thebibliography` | BibTeX + `oup-plain` |

**The decisive evidence is "Modern Large".** That design does not exist in
`NAR.cls`, whose only options are `a4`, `center`, and `fleqn`. The
Modern/Traditional/Contemporary × Large/Medium/Mediumone/Small matrix is
exclusively an `oup-authoring-template` feature (source 6). The instruction in
source 7 can therefore only refer to the OUP template.

Treat `NAR.cls` as legacy unless NAR's current site serves it — which we cannot
check, since the "File types → LaTeX" subsection is inside the truncated region
(§12).

**If we ever do have to switch,** the cost is contained by the architecture:
`sections/` contains zero styling commands, so only `main.tex` (18 front-matter
macros) and the four table environments would need rewriting.

### Useful intelligence from the legacy template

Even though we do not use it, it shows NAR house conventions the current
guidance has not given us:

1. **Conflict of interest heading.** It writes
   `\subsubsection{Conflict of interest statement.} None declared.` — NAR's own
   wording is "Conflict of interest statement", not "Conflict of interest".
   Ours should match. \TBD
2. **A CONCLUSION section exists** in its structure, between Discussion and
   Acknowledgements. Our manuscript has none. Confirm whether NAR expects one.
3. **Reference format**, visible in its `thebibliography`:
   `Author,A.B. and Author,C. (1992) Article title. *Abbreviated Journal Name*,
   **5**, 300--330.` — no space after the comma in initials, year in
   parentheses, journal abbreviated and italic, volume bold, en-dashed pages.
4. It hand-writes the bibliography **in citation order**, sidestepping the
   sorting problem in §11 entirely.
5. It uses `\footnote` for the corresponding author — supporting the reading in
   §13 that the "no footnotes" rule targets *content* footnotes only.

### Fonts are not set in the LaTeX

The class loads **no font package** — no `fontenc`, no Times/Helvetica/STIX
(verified: 32 `\RequirePackage` calls, none typographic). Compiled output
therefore renders in LaTeX's default Computer Modern and **will not look like a
published NAR paper**. This is expected: the author template fixes metrics,
layout, and structure, and OUP substitutes production fonts at typesetting.

Consequence for editing: text-mode Unicode is unreliable with no `fontenc`
declared, so Greek and thermodynamic symbols go through the math-mode macros
`\dfG`, `\drG`, `\drGo` defined in `main.tex` rather than being pasted in
literally.

The Arial / ≥7pt / 0.25–1pt rules in §3 are **artwork** requirements for the
image files. Nothing in the LaTeX enforces them.

---

## 8. Database requirements (not manuscript formatting)

- Must be **freely accessible without login**. Narrow exemptions for legally
  protected human data or acute funding constraints, agreed with the editor
  beforehand.
- HTTPS encouraged; required where sensitive data is handled.
- **URL persistence of at least five years** post-publication expected.
- Mobile and tablet accessibility encouraged, and "worth mentioning in the
  manuscript."
- Data availability must address "the formats and terms for data download."
- Supplementary material encouraged, must be complete at submission.

Source: [Database Issue Guidelines](https://academic.oup.com/nar/pages/Ms_Prep_Database).

---

## 9. Submission logistics

- **At least six** potential reviewers — "Authors must suggest at least six
  potential reviewers at submission." Add one extra for each reviewer you ask
  to exclude (exclusions go in the cover letter with a brief reason).
- Names, institutes, and email addresses.
  Independent (not recent collaborators), not from the same institution or city
  as any author. Omitting them "may delay handling."
- Disclosure of related recent or concurrent submissions required; duplicate
  submission "triggers automatic rejection and further sanctions."
- **Pre-submission technical checks:** NAR has partnered with Cactus for the
  Paperpal Preflight tool. Worth running before submission.
- Peer review is single-anonymized; manuscripts typically go to two reviewers.
- **Timing:** pre-submission enquiries to the Executive Editor by **1 July**;
  new-database manuscripts due **15 August**; **updates due 15 September**;
  nothing before 1 June.

Source: [Database Issue Guidelines](https://academic.oup.com/nar/pages/Ms_Prep_Database).

---

## 10. Initial submission

> "For the initial submission, we encourage you to submit a single .pdf file
> which includes the main text, references, tables, and figures. All figures and
> tables should be embedded in the text to facilitate reviewing."
>
> "Please upload supplementary data as separate file(s)."

**Do**

> - "Number all pages."
> - "Use embedded TrueType fonts in your Word document."
> - "Insert special characters using the Symbol font."
> - "Use single-column and single-spaced text (unless using LaTeX)"
> - "Submit a Graphical Abstract"

**Don't**

> - "Use line-numbering."
> - "Use footnotes."

Source 8 (see [Sources](#sources)) — NAR author guidelines, Manuscript
Preparation section, relayed by the corresponding author 2026-08-26. That
section could not be retrieved directly (see the retrieval notes).

### What this means for this project

| Requirement | Status |
|---|---|
| Single PDF, figures/tables embedded | **Satisfied by construction** — `latexmk -pdf main.tex` produces exactly this |
| Number all pages | **Satisfied** — the class puts `\thepage` in the running heads |
| Embedded TrueType fonts / Symbol font | **N/A** — both are Word instructions |
| Single-column, single-spaced | **Explicitly waived for LaTeX.** Keep the two-column Modern Large output; do not uncomment `\onecolumn` for the real submission |
| No line-numbering | **Satisfied** — `lineno` is not loaded (verified) |
| No footnotes | **Satisfied** — no `\footnote` anywhere in `sections/` or `main.tex` (verified). Keep it that way; use parenthetical text instead |
| Supplementary data as separate files | Nothing supplementary exists yet; the direction-sensitivity matrix is planned as supplementary |
| Graphical abstract | **Missing** — see §4 |

---

## 11. Known conflict: bibliography ordering

**This is unresolved and needs an answer from the editor or production.**

NAR requires (§6):

> "Cited in text by sequential number in order of appearance"

OUP's template manual instructs, for numbered style:

> "numbered citation style = `\bibliographystyle{oup-plain}`"

But `oup-plain.bst` **sorts the bibliography alphabetically**. Its `presort`
function builds a sort key from author, year, and title, then executes `SORT`
(line 1051 of `oup-plain.bst`):

```bibtex
FUNCTION {presort}
{ ... 'author.sort ... year field.or.null sortify ... title field.or.null ... }
ITERATE {presort}
SORT
```

With natbib in numeric mode and an alphabetically-sorted `.bst`, reference
numbers are assigned in alphabetical order. In-text citations therefore will
**not** ascend — you get (3), (1), (4) rather than (1), (2), (3).

This is not a misconfiguration on our side: `oup-plain` is exactly what OUP
documents for numbered style. The conflict is between OUP's supplied style file
and NAR's stated rule.

**Options, in order of preference:**

1. Ask the Executive Editor or production whether they renumber at typesetting.
   Most likely answer, and costs nothing to confirm.
2. Substitute an `unsrt`-derived `.bst` that preserves citation order. Deviates
   from OUP's documented instruction.
3. Use **`nar.bst`** (CTAN `/biblio/bibtex/contrib/misc/nar.bst`, v3.19,
   2011-12-23, "BibTeX style for Nucleic Acid Research"). Its header records
   that it "was created by Tom Schneider from unsrt.bst", and it contains
   **zero `SORT` commands** — so it preserves citation order, which is exactly
   what NAR asks for. Its output format also matches the reference style in the
   legacy NAR template (§7). Caveats: third-party rather than OUP-official,
   last updated 2011, and untested in combination with
   `oup-authoring-template.cls`. Verify a build before adopting.
4. Hand-author `\begin{thebibliography}` in citation order. The manual permits
   this ("The basic bibliography environment is accepted") but loses BibTeX.

Do not silently switch styles — record the decision here first.

---

## 12. Unverified — do not treat as settled

NAR's Author Guidelines "Manuscript preparation" section is **unreachable by
automated means**. Three URL forms were tried (plain, with the
`#section-13-7-10` anchor, and via a text-extraction proxy); direct `curl`
returns HTTP 403, proxies hit a CAPTCHA, and automated fetches truncate before
the section. Everything the corresponding author has relayed from it (§7, §10)
is marked as such.

Source 9 supplied most of the page but was **truncated at 50,000 characters,
before the "Preparing your manuscript" section**. From that page's own contents
list, the following remain unseen and are the highest-value gap:

> File types · LaTeX · Initial submission · Revision · **Article structure**
> (Title, Authors and affiliations, Graphical abstracts, Text Abstract,
> Introduction, Materials and Methods, Nomenclature conventions, Computer
> programs, Statistical analyses, Results, Discussion, Data Availability,
> Supplementary Data statement, Acknowledgements, Author Contributions
> Statement, Funding, Conflict of interest disclosure, **References**,
> **Tables**, **Figures**, BioRender, **Figure accessibility and alt text**,
> Supplementary Data)

**References, Tables and Figures in that list would very likely settle §11 and
the figure questions in §3.** Getting the rest of that page is the single most
useful remaining action.

The following are **guesses inherited from the OUP sample template or inferred
from the source manuscript**, not sourced requirements:

| # | Item | Where | Risk |
|---|---|---|---|
| 1 | `unnumsec` (unnumbered section heads) | `main.tex` class options | Copied from the sample's default line. No source says NAR wants unnumbered headings |
| 2 | `webpdf` | `main.tex` class options | Same. Manual defines it as "cropped paper size in the PDF output" |
| 3 | Section order and heading names | `main.tex` `\input` order | Follows `MANUSCRIPT.md`, not a NAR-specified order |
| 4 | Six keywords | `main.tex` `\keywords` | Invented. Unknown whether the Database issue uses keywords, or how many |
| 5 | `\appnotes{Database Issue}`, `\vol{00}`, `\issue{0}`, `\firstpage{1}` | `main.tex` metadata | Placeholder values from the sample |
| 6 | Back-matter **order** | `sections/author_contributions.tex` | The *set* is now confirmed (§13): Acknowledgements, Author Contributions, Funding and Conflict of interest disclosure all appear in NAR's Article structure list. Their required order and exact headings are still unseen |
| 7 | Abstract length | `sections/abstract.tex` | No limit found anywhere. Current draft is ~150 words |
| 8 | `table*` for the five-column tables | `sections/results_*.tex` | Layout judgement, not a requirement |
| 9 | Data Availability as a plain `\section` | `sections/data_availability.tex` | Heading and placement unconfirmed |
| 10 | NAR-specific figure overrides | §3 | The OUP artwork PDF says to use it "alongside any specific instructions provide on the website of the journal" — those instructions were never seen |

### Questions for the Executive Editor

Bundle these into one message; several are cheap to answer and unblock real work.

1. Does production renumber references, or should we supply a citation-order
   `.bst`? (§11)
2. Are there NAR-specific figure requirements that override the general OUP
   artwork guidance? (§3, item 10 above)
3. Should section headings be numbered or unnumbered? (item 1)
4. Is there an abstract word limit? (item 7)
5. Does the Database issue use keywords, and how many? (item 4)
6. The draft currently exceeds 4–6 typeset pages — is a longer update paper
   acceptable, or should we cut to fit? (§1)
7. Are there colour or page charges that scale with figure count? (The Charges
   section lists "Colour charges" and "Page charges" but the supplied text was
   truncated before their content.)

*Answered since first draft:* a conflict-of-interest disclosure **is** required
(§13), so that question has been removed from this list.

---

## 13. Authorship, ORCID, and mandatory disclosures

All quotes from source 9 (NAR Author Guidelines).

### ORCID is required

> "Submitting authors are required to provide an ORCID iD (Open Researcher and
> Contributor iD) at submission."

The OUP class supports this inline: `\author[1]{Name\ORCID{0000-0000-0000-0000}}`.

### Author Contributions statement is mandatory

> "The inclusion of an Author Contributions statement is mandatory for all
> articles, preferably with the original submission, but no later than at
> revision."

CRediT is encouraged. Non-author contributors go under Acknowledgements with
their contributions described.

### Conflict of interest disclosure is required

> "The Journal requires all authors to disclose any potential conflict of
> interest at the point of submission. It is the responsibility of the
> corresponding author to ensure that conflicts of interest of all authors are
> declared to the Journal."

"Conflict of interest disclosure" also appears as its own entry in the Article
structure contents list. This **resolves an earlier open question** — the
section is required, and the placeholder in
`sections/author_contributions.tex` is correct.

### AI disclosure — applies to this manuscript

> "The use of AI (for example, to help generate content or images, write code,
> process data, or for translation) should be disclosed both in cover letters to
> editors and in the Methods or Acknowledgements section of manuscripts."

> "Natural language processing tools driven by artificial intelligence (AI) do
> not qualify as authors, and the Journal will screen for them in author lists."

**This is a live obligation, not a hypothetical.** The LaTeX port of this
manuscript from `MANUSCRIPT.md`, and parts of the surrounding analysis tooling,
were produced with AI assistance. Decide with the co-authors what is disclosable
and write it into the cover letter and the Methods or Acknowledgements section.
Do not leave this to submission day.

### Equal contribution — and a footnote caveat

> "The joint authors should be identified by a dagger symbol and a footnote
> containing the statement '† Name1, Name2 and Name3 contributed equally to this
> work' should be added. The relative contributions of ALL authors must appear
> under Acknowledgements."

Note the tension with §10's "Don't use footnotes". The prohibition evidently
targets *content* footnotes in the body; the author-note apparatus (equal
contribution daggers, corresponding-author asterisks) is explicitly required and
is handled by the class's `\corresp` and author-mark machinery, not by
`\footnote`. **Keep `\footnote` out of the body regardless.**

---

## 14. Computational resources — the database itself

This is a Database paper, so NAR's computational-resource policy is a
requirement on the resource, not just the manuscript. All quotes from source 9.

The journal requires that any published database, webserver, or dataset be:

> "freely available over the internet, without login or registration, and"
> "updated or at least maintained in a fully functional form, ideally at the same
> URL for at least 5 years."

Databases specifically:

> "If any part of the database (e.g. the one that deals with the user-submitted
> data) needs to be password-protected, only the freely available part will be
> considered by the reviewers."

> "Authors are encouraged, but not required, to make the contents of their
> databases freely available as flat or relational files upon request."

Software:

> "Authors must ensure that the software is available for a full 5 years
> following publication, preferably through a download link on a stable URL or in
> a public code repository such as GitHub."

### Code deposit — action required at revision

> "you will be required to deposit the code to Zenodo or FigShare and provide the
> permanent DOI in the Data Availability statement at revision or acceptance.
> These repositories will retain the version of the code that was used in your
> paper."

A GitHub URL alone is **not** sufficient. `sections/data_availability.tex`
currently cites the repository and a planned PyPi package; it needs a Zenodo or
FigShare DOI by revision.

Data availability statements are also mandatory:

> "The inclusion of a data availability statement is a requirement for papers
> published in the Journal."


---

## Sources

| # | Source | URL | Retrieved |
|---|---|---|---|
| 1 | NAR Database Issue Guidelines | https://academic.oup.com/nar/pages/Ms_Prep_Database | 2026-08-26 |
| 2 | OUP — *Guidance for preparing artwork* (PDF) | https://static.primary.prod.gcms.the-infra.com/static/site/journals/document/images-author-guidance.pdf?node=1bf05d0b2fbd9c529a23&version=490455:30c2211aa70bba63a5ee | 2026-08-26 |
| 3 | NAR Web Server Issue Guidelines (cross-check on length) | https://academic.oup.com/nar/pages/Submission_Webserver | 2026-08-26 |
| 4 | NAR Author Guidelines | https://academic.oup.com/nar/pages/author-guidelines | 2026-08-26 |
| 5 | NAR Methods Guidelines | https://academic.oup.com/nar/pages/methods-guidelines | 2026-08-26 |
| 6 | `oup-authoring-template` v1.5 on CTAN | https://ctan.org/pkg/oup-authoring-template | 2026-08-26 |
| 7 | NAR LaTeX design requirement ('Modern Large') — relayed by the corresponding author; originates from OUP's *Preparing and submitting your manuscript* page, which returned navigation-only content on automated fetch | https://academic.oup.com/journals/pages/authors/preparing_your_manuscript | 2026-08-26 |
| 8 | NAR author guidelines, Manuscript Preparation section (initial-submission rules, Do/Don't list) — relayed by the corresponding author; not retrievable by automated means | https://academic.oup.com/nar/pages/author-guidelines#section-13-7-10 | 2026-08-26 |
| 9 | NAR Author Guidelines — full page text (scope, peer review, ethics, authorship, ORCID, AI disclosure, data/code availability, computational resources, charges), supplied by the corresponding author. **Truncated at 50,000 characters before the "Preparing your manuscript" section**, so Article structure / References / Tables / Figures remain unseen | https://academic.oup.com/nar/pages/author-guidelines | 2026-08-26 |
| 10 | `nar.bst` v3.19 — BibTeX style for Nucleic Acid Research (candidate fix for §11) | https://ctan.org/pkg/nar | 2026-08-26 |
| 11 | OUP General Template on Overleaf (same class as source 6) | https://www.overleaf.com/latex/templates/oup-general-template/ybpypwncdxyb | 2026-08-26 |

A local copy of source 2 is **not** committed here — it is OUP copyright. Re-download
from the URL above if the quotes need re-checking.

### Retrieval notes

Sources 4 and 5 were reached but returned **no** length, figure-count, or
figure-format content — the Author Guidelines' "Manuscript preparation" section
sits behind Cloudflare (direct `curl` returns HTTP 403, reader proxies hit the
CAPTCHA, and automated fetches truncate before reaching it). Three URL forms
were tried, including the `#section-13-7-10` anchor; all three truncate at the
same point. Sources 7 and 8 come from that section and were relayed by the
corresponding author. **Anything else in that section remains unknown** — see
§12. The figure
specification in §3 therefore comes from OUP's central artwork PDF (source 2),
not from a NAR-specific page. Per that PDF's own caveat, **check the NAR site for
journal-specific overrides** before final artwork submission.

---

## Open items for this manuscript

| # | Item | Where | Blocked on |
|---|---|---|---|
| 1 | Title should begin with the database name, not "The" | `main.tex` title block | author decision |
| 2 | Graphical abstract does not exist yet | separate file at submission | design work |
| 3 | Draft is over the 4–6 page budget — compile and measure | whole document | a TeX build |
| 4 | `[VERIFY]` notes in the bibliography unresolved | `references.bib` | lookup |
| 5 | Data availability lacks download formats, 5-year persistence, mobile note | `sections/data_availability.tex` | drafting |
| 6 | Six referees to nominate | ScholarOne, at submission | author decision |
| 7 | Figures as vector PDF/EPS, Arial, ≥7pt, 0.25–1pt lines | figure generator scripts | figure work |
| 8 | Bibliography ordering conflict | `main.tex` `\bibliographystyle` | **editor** (§11) |
| 9 | Ten unverified formatting assumptions | various | **editor** (§12) |
| 10 | Keep the document free of footnotes and line numbers | all `sections/` | ongoing discipline |

## Compliance already verified

Checked against the source, not assumed:

- No `\footnote` anywhere in `main.tex` or `sections/` — NAR prohibits footnotes
- `lineno` not loaded — NAR prohibits line numbering
- Pages numbered via the class running heads
- Build produces a single PDF with tables embedded
- All four tables' column counts match their specs
- No bare underscores across 84 `\file`/`\texttt`/`\cpd` uses
- No non-ASCII characters left from the markdown conversion
- All 25 `\input` targets resolve; braces balanced; all cite keys and refs resolve
