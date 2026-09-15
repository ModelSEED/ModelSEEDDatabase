# Thermodynamics pipeline: dependencies that are not redistributed

Two inputs the thermodynamics pipeline consumes are **not shipped in this
repository**. Both are reachable by citation, but neither can be regenerated
from what is released here. This file exists so that a reader reconciling the
released tables against their deposited sources knows why a referenced input is
absent.

## 1. Digitized measured pKa values (IUPAC)

Zheng, J. and Lafontant-Joseph, O. *IUPAC/Dissociation-Constants: v2.3b*,
Zenodo, 2025. <https://doi.org/10.5281/zenodo.15375522>
Digitized from Serjeant & Dempsey (1979) and Perrin (1965, 1972 suppl.).

Released under **CC-BY-NC-4.0**, which is incompatible with this repository's
MIT licence, so the file is consumed by the pipeline but deliberately not
committed (it is listed in `.gitignore`). A reader reconciling released tables
against the deposited values will therefore find entries sourced from a file
that is not present. Download it from the Zenodo DOI above to reproduce that
layer.

## 2. eQuilibrator compound cache

The cache is a **pinned public release** rather than an artefact rebuilt here.
This makes the protonation layer reproducible by citation but not regenerable
without a licensed tool.

---

Everything else needed to regenerate published values — including the
eQuilibrator cache rebuild scripts — is under `Scripts/Thermodynamics/`.

*Moved out of the manuscript's Data Availability section, 2026-09-15.*
