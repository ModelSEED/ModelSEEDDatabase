# eQuilibrator regeneration pipeline

The chain that rebuilds the eQuilibrator compound cache from ModelSEED
structures, refits component contribution against it, and writes the energy
tables that `../Rerun_Thermodynamics.sh` then ingests.

**Two stages, and only the second is in `Rerun_Thermodynamics.sh`:**

| stage | what it does | needs |
|---|---|---|
| this directory | cache → refit → `ModelSEED_{Reaction,Compound}_Energies.tsv` | the eQuilibrator working tree, ~4 GB, and hours |
| `../Rerun_Thermodynamics.sh` | those committed tables → the reaction/compound JSON | this repository only |

Most people want the second. Run this one only when the structures, the pKa
layer or the training data have changed.

## Why it is here

It was unversioned. The scripts lived in one author's `eQuilibrator/tools`
directory, and two already-committed scripts imported from it, so the
dependency was live rather than hypothetical. A release that cannot be
regenerated from the repository is not reproducible, and several working
corrections existed only as loose files on a single filesystem.

## What it needs

    EQUILIBRATOR_DIR   the eQuilibrator working tree (default: the author's
                       path, which is almost certainly wrong for you)

That tree holds the compound caches, the fitted parameter files and the pinned
`msd_dev` structure export. None of it is committed here — the caches alone are
gigabytes. Paths *into this repository* are derived from the script location and
need no configuration.

    export EQUILIBRATOR_DIR=/path/to/eQuilibrator

## Order

```
build_seed_mapping.py         seed:cpd → structure → cache-compound mapping
must_carry_over.py            compounds that genuinely cannot be created
write_resolved_pkas.py        the pKa cascade, materialised as one table
build_modelseed_cache.py      the cache, built from ModelSEED structures
train_path_b.py               refit component contribution against that cache
generate_modelseed_energies.py            reaction ΔG′°
generate_modelseed_compound_energies.py   compound ΔfG′°
update_modelseed_thermodynamics.py        write both into the JSON
verify_regen.py               diff a regeneration against what is committed
```

`modelseed_pkas.py` is the shared library the first four import — structures,
the pKa sources and the resolution order.

## The pKa resolution order

`alberty → iupac → cache → marvin → molgpka`, with the two ChemAxon-derived
tiers **gated**: they are consulted only where MolGpKa provably cannot answer,
detected mechanically as a repeated predicted value, which is the signature of
two chemically equivalent sites returning the same number. That holds for 6,664
of 32,201 compounds.

Gating matters beyond provenance. Substituting a compound's source changes its
site count above pH 7, and if a reaction partner keeps the old enumeration the
two no longer cancel — that is what moved `rxn08789` by 370 kcal/mol when the
cache tier was applied unconditionally. Every needless substitution is a chance
to break a reaction.

Two further guards in `modelseed_pkas.py`: a ladder the magnesium guard would
refuse falls through to the next tier instead of being installed and silently
skipped, and monatomic ions assert an empty list rather than letting a fallback
invent the ChemAxon pKa-3.09 artefact.

## Guards worth knowing about

**`update_modelseed_thermodynamics.py` refuses inputs that did not come from the
cache being shipped.** It compares the provenance header of each input table
against the shipped cache *by content*, not by path — a cache is routinely
promoted under a new name between generating a table and shipping it. On
2026-09-04 both inputs were five days stale, the run reported "23,800 updated"
and the entry count rose, and the JSON received pre-rebuild numbers with nothing
in the summary to indicate it. `--allow-stale` downgrades the refusal to a
warning.

**The pKa window is −2..16, not 0..14.** eQuilibrator counts sites above the
reported pH to place the major microspecies, so a site at pKa 14.9 is protonated
at pH 7 and carries a proton. `equilibrator-assets` clipped user-supplied pKas
to 0..14 regardless of what the caller asked for; the fix is
`EQUILIBRATOR_PKA_MIN` / `EQUILIBRATOR_PKA_MAX` in that package. Measured
effect on the energies at pH 7: **none** — the major species is anchored to the
compound's `atom_bag`, so adding states at the extremes does not move it. It
does make the stored ladders complete, which is why it is on.

## Known limits

* **133 reactions carry an unreliable protonation model**, listed in
  `Biochemistry/Thermodynamics/ProtonationEvidence/unbalanced_enumeration_reactions.tsv`.
  Large glycans whose partners could not move with them. Neither the current nor
  the previous value can be checked against experiment; none has a TECRDB match.
* **The rebuild is not measurably better against experiment.** All build
  variants score MAE 0.74 against the 802 stereo-exact TECRDB anchors, because
  none of the reactions the rebuild changes has an experimental measurement. The
  case for it is correctness of construction, not demonstrated accuracy.
