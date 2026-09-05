"""Path B: retrain Component Contribution against our own compound cache.

The training reactions are stored as accession strings (kegg:C00002 + ...) and
resolved at train time through ``ccache.get_compound``. Trainer.train() then
indexes train_S and train_G by ``compound.id``. So the compound ids baked into
the resulting cc_params come from whichever cache is passed in -- which is what
makes a ModelSEED-anchored model possible without touching the algorithm.

Stage 1 (--cache path_a) reproduces upstream: the kegg:/metanetx: accessions
resolve exactly as they do in the shipped cache, so the parameters should match
what Zenodo ships. That is the control -- if it does not reproduce, nothing
downstream is trustworthy.
"""

import os
import argparse
import logging
import sys
import time
from pathlib import Path

# Relocated from the eQuilibrator working tree 2026-09-05. ROOT is that
# tree -- it holds the caches and fitted parameters, which are gigabytes and
# are not in this repository -- so it is named by environment variable rather
# than derived from this file's location.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

DEFAULT_CACHE = ROOT / "data" / "modelseed_cache" / "compounds.sqlite"
DEFAULT_OUT = ROOT / "data" / "cc_params_path_b.npz"



def fixed_tecrdb() -> Path:
    """Work around a silent pandas-3 breakage in ``read_tecrdb``.

    component-contribution 0.7.0 applies its documented defaults with

        tecr_df.ionic_strength.fillna(0.25, inplace=True)
        tecr_df.p_mg.fillna(14, inplace=True)

    Under pandas >= 3 copy-on-write those operate on a temporary Series and do
    nothing at all -- pandas only emits a ChainedAssignmentError *warning*. In
    this TECRDB snapshot 75% of rows have no ionic strength and 66% have no
    pMg, so three quarters of the training set would carry NaN conditions
    straight into the reverse transform.

    Rather than patch site-packages, write a copy with the defaults already
    applied and hand it to the factory explicitly.
    """
    import pandas as pd
    from equilibrator_cache.zenodo import get_cached_filepath
    from component_contribution import TRAINING_DATA_TECR_SETTINGS

    src = get_cached_filepath(TRAINING_DATA_TECR_SETTINGS)
    df = pd.read_csv(src)
    before = (int(df.ionic_strength.isna().sum()), int(df.p_mg.isna().sum()))
    df["ionic_strength"] = df["ionic_strength"].fillna(0.25)
    df["p_mg"] = df["p_mg"].fillna(14)
    out = ROOT / "data" / "TECRDB_defaults_applied.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"TECRDB defaults applied: ionic_strength {before[0]} NaN -> "
          f"{int(df.ionic_strength.isna().sum())}, "
          f"p_mg {before[1]} -> {int(df.p_mg.isna().sum())}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--tecrdb", type=Path, default=None,
                    help="use this TECRDB CSV instead of the Zenodo one")
    args = ap.parse_args()

    if args.quiet:
        logging.getLogger("component_contribution").setLevel(logging.CRITICAL)
        logging.getLogger("equilibrator_cache").setLevel(logging.CRITICAL)

    from component_contribution.training_data import FullTrainingDataFactory
    from component_contribution.trainer import Trainer
    from equilibrator_cache.api import create_compound_cache_from_sqlite_file

    print(f"cache: {args.cache}")
    ccache = create_compound_cache_from_sqlite_file(args.cache)

    tecrdb_path = args.tecrdb if args.tecrdb else fixed_tecrdb()

    print("building training data (parses TECRDB + formation + redox) ...")
    t0 = time.time()
    td = FullTrainingDataFactory(ccache=ccache).make(tecrdb_path=tecrdb_path)
    print(f"  built in {time.time() - t0:.0f}s")
    print(f"  reactions                 : {td.stoichiometric_matrix.shape[1]}")
    print(f"  compounds                 : {len(td.compounds)}")
    print(f"  decomposable              : {len(td.decomposable_compounds)}")
    print(f"  non-decomposable          : {len(td.non_decomposable_compounds)}")
    print(f"  group definitions         : {td.group_summary.shape[0]}")

    print("\ntraining ...")
    t0 = time.time()
    params = Trainer.train(td)
    print(f"  trained in {time.time() - t0:.0f}s")
    print(f"  dimensions:\n{params.dimensions}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    params.to_npz(args.out)
    print(f"\nwrote {args.out} ({args.out.stat().st_size / 1e6:.1f} MB)")

    ids = list(params.train_G.index)
    print(f"  train_G indexed by {len(ids)} compound ids, range {min(ids)}-{max(ids)}")


if __name__ == "__main__":
    main()
