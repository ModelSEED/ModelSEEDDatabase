# Thermodynamics pipeline: dependencies that are not redistributed

Three inputs the thermodynamics pipeline consumes are **not shipped in this
repository**. All are reachable by citation, but none can be regenerated
from what is released here. This file exists so that a reader reconciling the
released tables against their deposited sources knows why a referenced input is
absent.

## 1. Digitized measured pKa values (IUPAC)

Zheng, J. and Lafontant-Joseph, O. *IUPAC/Dissociation-Constants: v2.3b*,
Zenodo, 2025. <https://doi.org/10.5281/zenodo.15375522>
Digitized from Serjeant & Dempsey (1979) and Perrin (1965, 1972 suppl.).

Released under **CC-BY-NC-4.0**. The non-commercial clause is incompatible with
this repository's Creative Commons Attribution licence, so the file is consumed
by the pipeline but deliberately not committed (it is listed in `.gitignore`). A reader reconciling released tables
against the deposited values will therefore find entries sourced from a file
that is not present. Download it from the Zenodo DOI above to reproduce that
layer.

## 2. eQuilibrator compound cache

The cache is a **pinned public release** rather than an artefact rebuilt here.
This makes the protonation layer reproducible by citation but not regenerable
without a licensed tool.

## 3. dGbyG code and model weights

Fan, W. *et al.* Unraveling principles of thermodynamics for genome-scale
metabolic networks using graph neural networks. *Cell Systems* (2025).
<https://doi.org/10.1016/j.cels.2025.101393>. Code and weights:
<https://github.com/f-wc/dGbyG> (MIT), at the commit recorded in the header
line of `Biochemistry/Thermodynamics/dGbyG/ModelSEED_*_Energies.tsv`.

`Generate_dGbyG_Energies.py` runs that checkout's shipped ensemble
(`models/mpnn_A139_B23_E300_L2_v2`, 100 heads) over ModelSEED's structures to
produce the two staged tables. The checkout and its weights are not copied
here; the staged tables are, so the database can be rebuilt from them without
dGbyG, and regenerating them needs the checkout and its Python dependencies
(dGbyG's `environment.yml`).

---

Everything else needed to regenerate published values — including the
eQuilibrator cache rebuild scripts — is under `Scripts/Thermodynamics/`.

*Moved out of the manuscript's Data Availability section, 2026-09-15.*
