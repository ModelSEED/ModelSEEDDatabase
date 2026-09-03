"""Path resolution for the protonation-evidence scripts.

The repository root is derived from this file's location, so the in-tree paths
need no configuration. Everything outside the repository -- the eQuilibrator
working tree, the caches it builds, and the intermediate tables from the pKa
experiment -- is resolved from an environment variable with a documented
default, because those live on the analysis host rather than in version
control.

    EQUILIBRATOR_DIR    checkout of the eQuilibrator working tree
    EQUILIBRATOR_CACHE  the pristine upstream (Zenodo) compound cache
    MSD_WORKDIR         scratch directory holding the pKa-experiment tables
    THERMO_CACHE        the rebuilt ModelSEED cache being graded

Missing inputs raise on use rather than silently degrading: a grading run that
quietly skipped half its inputs would report a flattering distribution.
"""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]

EQ = Path(os.environ.get(
    "EQUILIBRATOR_DIR", "/scratch/seaver/Claude_Projects/eQuilibrator"))
WORK = Path(os.environ.get(
    "MSD_WORKDIR", "/scratch/seaver/Claude_Projects/MSD_Structures/pka_experiment"))
ZENODO_CACHE = Path(os.environ.get(
    "EQUILIBRATOR_CACHE", Path.home() / ".cache/equilibrator/compounds.sqlite"))
THERMO_CACHE = Path(os.environ.get(
    "THERMO_CACHE", EQ / "data/cache_final/compounds.sqlite"))

# in-repo, always present
STRUCTURES = REPO / "Biochemistry/Structures/Unique_ModelSEED_Structures.txt"
PKA_DIR = REPO / "Biochemistry/Structures/ModelSEED/pkas"
DATA_OUT = REPO / "Biochemistry/Thermodynamics/ProtonationEvidence"

# produced by the pKa experiment, outside the repo
RANKED = WORK / "polyprotic_ranked.tsv"
REACTIONS = WORK / "reactions_final.tsv"
RESOLVED = EQ / "data/resolved_pkas.tsv"
PRIORITY = REPO.parent / "priority/template_reactions.txt"

IUPAC_ZIP = Path(os.environ.get(
    "IUPAC_PKA_ZIP", Path.home() / ".cache/iupac_pka.zip"))


def require(path, what):
    """Fail loudly. A silently-skipped input produces a flattering result."""
    if not Path(path).exists():
        raise SystemExit(
            f"missing {what}: {path}\n"
            f"set the corresponding environment variable (see paths.py) "
            f"or generate the file first.")
    return Path(path)
