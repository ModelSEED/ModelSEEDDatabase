# KBase biochemistry object — frozen at the 2020 release

These scripts are **legacy** and are deliberately not maintained. The KBase
narrative is no longer updated with the biochemistry as a matter of course; it
is regenerated only on request.

Everything here is restored to its state at the 2020 publication
(Seaver *et al.*, NAR 2021;49(D1):D575–D588, `10.1093/nar/gkaa746`):

| file | restored from |
|---|---|
| `Build_Biochemistry_Object.py` | `b83c6fcd` (2020-08-24), the last change before publication |
| `Load_WS_Biochemistry.py`, `Load_MS_Biochemistry.pl` | unchanged since 2020 |

The input skeleton `Objects/Base_Biochemistry.json` (id `MSD_v1.0`, 16
compartments, 306 cues) was last committed at `43026f61` (2020-08-19) and is
already at its 2020 state. It stays in `Objects/` because that directory is
shared with the active disambiguation curation scripts, which write working
files there — renaming it would break them. The generated output
(`MSD_v1.0_Biochem.json`) has never been tracked.

## It will not run against the current database

This is intentional, and worth stating plainly rather than discovering:

- it reads `Biochemistry/compounds.json` and `Biochemistry/reactions.json`, the
  pre-split monolithic files, which were replaced by the sharded
  `compound_NN.json` / `reaction_NN.json` layout;
- it reads the flat `deltag`, `deltagerr`, `pka` and `pkb` fields, which were
  removed on 2026-09-11 — energies now live per source under `thermodynamics`
  and protonation per tool under `pkas`.

To regenerate a KBase object from current data, the mapping has to be rewritten
against those two structures. A version that did so briefly existed on
2026-09-11 and was reverted in favour of this freeze; see the git history of
`Scripts/Release/Build_Biochemistry_Object.py` if it is ever wanted back.
