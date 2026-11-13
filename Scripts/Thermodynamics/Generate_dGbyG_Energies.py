#!/usr/bin/env python
"""Regenerate the dGbyG energy tables from ModelSEED's own structures.

dGbyG (Fan et al. 2025, Cell Systems, doi:10.1016/j.cels.2025.101393;
https://github.com/f-wc/dGbyG) is a graph-neural-network ensemble that predicts
standard transformed formation energies. This is the UPSTREAM stage for that
source: it runs the shipped dGbyG ensemble over every compound with a SMILES
and sums the result over every reaction. Update_{Compound,Reaction}_dGbyG_
Energies.py then install the tables into the JSON, exactly as the eQuilibrator
updaters install theirs.

  input   Biochemistry/compound_*.json  (smiles)
          Biochemistry/reaction_*.json  (stoichiometry)
          a dGbyG checkout (--dgbyg), its shipped ensemble
          models/mpnn_A139_B23_E300_L2_v2 (100 heads)
  output  Biochemistry/Thermodynamics/dGbyG/ModelSEED_Compound_Energies.tsv
          Biochemistry/Thermodynamics/dGbyG/ModelSEED_Reaction_Energies.tsv

One row per compound / reaction record, live and obsolete, with a `status`
column saying why a record has no energy; only `ok` rows carry one. The first
line records the dGbyG commit, the model actually loaded, its head count and
the conditions, and must not be dropped from the file.

CONDITIONS. dGbyG's native condition: pH 7.0, ionic strength 0.25 M, pMg 14,
298.15 K -- the condition its training data were transformed to. Not pMg 3.0
like the eQuilibrator tables: dGbyG's Legendre transform has no Mg-binding
constants (it moves only explicit Mg atoms), so re-transforming would not model
Mg binding anyway. No pKa values are needed at the native condition.

COMPOUNDS. dGbyG's own normalisation (rdMolStandardize Normalize + Uncharger)
of each ModelSEED SMILES; the value is the mean over the 100 heads and the
uncertainty their standard deviation. Multi-fragment SMILES (salts,
counter-ions) are predicted as given. Statuses:
  no structure        no SMILES in the compound record
  R group             the SMILES carries an attachment point ('*'). dGbyG would
                      featurise the dummy atom and return a number, but it is an
                      extrapolation outside anything the model was trained on.
  not parsed          dGbyG could not build or normalise a molecule from the
                      SMILES (no molecule, or RDKit raised, e.g. a valence error)
  prediction failed   dGbyG raised while featurising or predicting

REACTIONS. Stoichiometry is netted by compound across compartments -- what the
eQuilibrator and dGPredictor tables also score -- and H+ (cpd00067) is dropped:
its transformed formation energy is 0 at fixed pH, and dGbyG's normalisation
changes the protonation state of many species from ModelSEED's pH-7 form, so
ModelSEED's H+ count no longer balances what dGbyG sees. The value is the mean
over heads of the SUMMED per-head vector, and the uncertainty is the standard
deviation of that summed vector -- dGbyG's own Reaction.standard_dGr_prime.
Compound uncertainties are NOT added in quadrature: the heads' errors are
strongly correlated (ATP's compound sigma is 5.7 kcal/mol, ATP hydrolysis's
0.8). Statuses, in order of precedence:
  no stoichiometry     the record has none
  translocation only   every participant cancels across compartments
  missing structure    a participant has no SMILES
  compound failed      a participant is 'not parsed' or 'prediction failed'
  R group              a participant carries an R group (see above)
  unbalanced           fails dGbyG's balance check (Reaction.is_balanced:
                       ignore_H, ignore_H_ion) on dGbyG-normalised structures.
                       Heavy atoms and water are enforced; H is not, and a
                       charge difference that matches an H difference is
                       absorbed as H+, so electron-unbalanced half-reactions
                       pass. dGbyG's Reaction.balance() is meant to add water
                       for an oxygen deficit but never does (src/dGbyG/api.py
                       tests is_balanced on the unmodified reaction), so dGbyG
                       itself returns NaN for these.
  degenerate           every head gives the same value (sigma < 0.005 kcal/mol,
                       i.e. 0.00 at the stored precision): the two sides are the
                       same graph up to stereochemistry -- dGbyG has no chirality
                       or E/Z feature -- or the reaction leaves every atom's
                       2-bond neighbourhood unchanged (a 2-layer message-passing
                       network cannot see past that; e.g. NAD + NADPH = NADH +
                       NADP). The value is exactly 0 and its zero uncertainty is
                       a blind spot, not certainty, so no energy is shipped.

PRECISION. kJ/mol to 3 decimals (as the eQuilibrator tables); kcal/mol rounded
to the 2 decimals the database stores, from the unrounded kJ value, so the
updaters' 2-decimal formatting is the identity and the operator they compute at
write time sees exactly the stored numbers.

Inference runs on CPU in --workers spawned processes (each loads the ensemble
once). Keep it on CPU: GPU float32 arithmetic can move a value across a
2-decimal rounding boundary, and the committed tables were produced on CPU.

REQUIRES. A dGbyG checkout (https://github.com/f-wc/dGbyG) with its shipped
model weights, and dGbyG's Python dependencies (its environment.yml): torch,
torch_geometric, rdkit, numpy, pandas, JPype1, pubchempy, biopython, tqdm,
requests and portalocker.
"""

if __name__ == "__main__":
    # Argument guard -- see "The argument guard" in Scripts/README.md.
    import argparse as _argparse
    _parser = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter)
    _parser.add_argument(
        "--dgbyg", required=True,
        help="path to a dGbyG checkout (the directory holding src/ and models/)")
    _parser.add_argument(
        "--workers", type=int, default=8,
        help="inference processes (default 8); each loads the 100-head ensemble")
    _parser.add_argument(
        "--threads", type=int, default=2,
        help="torch threads per process (default 2)")
    _ARGS = _parser.parse_args()


import csv
import glob
import json
import multiprocessing as mp
import os
import subprocess
import sys
import warnings
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIOCHEM = os.path.normpath(os.path.join(HERE, "..", "..", "Biochemistry"))
OUT_DIR = os.path.join(BIOCHEM, "Thermodynamics", "dGbyG")
KJ_PER_KCAL = 4.184
PROTON = "cpd00067"
MODEL = "mpnn_A139_B23_E300_L2_v2"
CONDITION = "p_h=7.0 ionic_strength=0.25M p_mg=14 T=298.15K"
DEGENERATE_SIGMA = 0.005     # kcal/mol: rounds to 0.00 at the stored precision

_worker = {}


def _quiet():
    warnings.filterwarnings("ignore")
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")


def _init_worker(dgbyg_dir, threads):
    """Load dGbyG and its ensemble once per process, on CPU."""
    _quiet()
    import torch
    torch.set_num_threads(threads)
    sys.path.insert(0, os.path.join(dgbyg_dir, "src"))
    import dGbyG.api as api
    from dGbyG.model.inference import Inference_Model
    api.model_cache[api.infer_model_path] = Inference_Model(api.infer_model_path, device="cpu")
    _worker["api"] = api


def _predict(item):
    """(cpd, smiles) -> (cpd, status, atom bag, per-head kJ/mol or None)."""
    cpd, smiles = item
    api = _worker["api"]
    try:
        comp = api.Compound(smiles, "smiles")
    except Exception:
        return cpd, "not parsed", None, None
    if comp.mol is None:
        return cpd, "not parsed", None, None
    try:
        bag = comp.atom_bag
        heads = np.asarray(comp.standard_dGf_prime_list, dtype=np.float32)
    except Exception:
        return cpd, "prediction failed", None, None
    if heads.ndim != 1 or not np.isfinite(heads).all():
        return cpd, "prediction failed", None, None
    return cpd, "ok", bag, heads


def preflight(dgbyg_dir):
    """Fail fast, in this process, on anything that would otherwise make every
    spawned worker die in its initializer. Returns (model name, head count)
    as dGbyG itself resolves them, for the provenance header."""
    if not os.path.isdir(os.path.join(dgbyg_dir, "src", "dGbyG")):
        sys.exit(f"--dgbyg {dgbyg_dir}: no src/dGbyG package there")
    _quiet()
    sys.path.insert(0, os.path.join(dgbyg_dir, "src"))
    try:
        import dGbyG.api as api
    except ImportError as e:
        sys.exit(f"cannot import dGbyG from {dgbyg_dir} ({e}); see REQUIRES in --help")
    paths = api.infer_model_path
    if len(paths) != 1 or os.path.basename(os.path.normpath(paths[0])) != MODEL:
        sys.exit(f"dGbyG resolves its model to {paths}, not {MODEL}; update MODEL "
                 f"and the README if that is intended")
    heads = sorted(glob.glob(os.path.join(paths[0], "*.pt")))
    if len(heads) != len(os.listdir(paths[0])):
        sys.exit(f"{paths[0]} holds files other than *.pt; dGbyG would load them as heads")
    return MODEL, len(heads)


def is_balanced(net, bags):
    """dGbyG's reaction_utils.is_balanced on atom bags, with the flags its
    Reaction class uses (ignore_H, ignore_H_ion)."""
    diff = defaultdict(float)
    for cpd, coeff in net.items():
        for atom, num in bags[cpd].items():
            diff[atom] += num * coeff
    if diff.get("charge", 0) * diff.get("H", 0) <= 0:
        h_ion = 0
    elif diff["charge"] < 0:
        h_ion = -max(diff["charge"], diff["H"])
    else:
        h_ion = -min(diff["charge"], diff["H"])
    diff["charge"] = diff.get("charge", 0) + h_ion
    diff["H"] = 0
    return not any(abs(v) > 1e-9 for v in diff.values())


def dgbyg_commit(dgbyg_dir):
    """HEAD of the checkout, marked '+modified' if anything dGbyG loads -- its
    source or any file in models/, tracked or not -- differs from it."""
    try:
        sha = subprocess.run(["git", "-C", dgbyg_dir, "rev-parse", "HEAD"], capture_output=True,
                             text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", dgbyg_dir, "status", "--porcelain", "--untracked-files=all",
                                "--", "src", "models"], capture_output=True, text=True, check=True).stdout
        return sha + ("+modified" if dirty.strip() else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def load(kind):
    records = []
    for shard in sorted(glob.glob(os.path.join(BIOCHEM, f"{kind}_[0-9][0-9].json"))):
        with open(shard) as fh:
            records.extend(json.load(fh))
    return records


def kcal2(kj):
    """kJ/mol -> kcal/mol, rounded to the 2 decimals the database stores."""
    return f"{round(kj / KJ_PER_KCAL, 2):.2f}"


def main(args):
    dgbyg_dir = os.path.abspath(args.dgbyg)
    model, n_heads = preflight(dgbyg_dir)
    header = (f"# dgbyg=https://github.com/f-wc/dGbyG@{dgbyg_commit(dgbyg_dir)} model={model} "
              f"heads={n_heads} {CONDITION}\n")
    compounds, reactions = load("compound"), load("reaction")
    print(f"{len(compounds):,} compound records, {len(reactions):,} reaction records")

    status, bags, heads = {}, {}, {}
    work = []
    for c in compounds:
        smiles = c.get("smiles") or ""
        if not smiles:
            status[c["id"]] = "no structure"
        elif "*" in smiles:
            status[c["id"]] = "R group"
        else:
            work.append((c["id"], smiles))
    print(f"predicting {len(work):,} compounds on {args.workers} CPU workers", flush=True)
    # ProcessPoolExecutor, not multiprocessing.Pool: a worker that dies in its
    # initializer breaks the executor (and this run) instead of being respawned
    # forever.
    with ProcessPoolExecutor(args.workers, mp_context=mp.get_context("spawn"),
                             initializer=_init_worker, initargs=(dgbyg_dir, args.threads)) as pool:
        for n, (cpd, st, bag, h) in enumerate(pool.map(_predict, work, chunksize=32), 1):
            status[cpd] = st
            if h is not None:
                if h.shape[0] != n_heads:
                    sys.exit(f"{cpd}: {h.shape[0]} head values, expected {n_heads}")
                bags[cpd], heads[cpd] = bag, h
            if n % 5000 == 0:
                print(f"  {n:,}/{len(work):,}", flush=True)

    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "ModelSEED_Compound_Energies.tsv")
    with open(path + ".tmp", "w", newline="") as fh:
        fh.write(header)
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["compound_id", "name", "formula", "status", "dgf_prime_kJ_per_mol",
                    "uncertainty_kJ_per_mol", "dgf_prime_kcal_per_mol", "uncertainty_kcal_per_mol"])
        for c in sorted(compounds, key=lambda r: r["id"]):
            st = status[c["id"]]
            vals = [""] * 4
            if st == "ok":
                h = heads[c["id"]].astype(np.float64)
                m, s = h.mean(), h.std()
                vals = [f"{m:.3f}", f"{s:.3f}", kcal2(m), kcal2(s)]
            w.writerow([c["id"], c.get("name", ""), c.get("formula", ""), st] + vals)
    os.replace(path + ".tmp", path)
    print("compounds:", dict(Counter(status.values())))

    rstatus = Counter()
    path = os.path.join(OUT_DIR, "ModelSEED_Reaction_Energies.tsv")
    with open(path + ".tmp", "w", newline="") as fh:
        fh.write(header)
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["reaction_id", "name", "status", "dg_prime_kJ_per_mol", "uncertainty_kJ_per_mol",
                    "dg_prime_kcal_per_mol", "uncertainty_kcal_per_mol"])
        for rxn in sorted(reactions, key=lambda r: r["id"]):
            stoich = rxn.get("stoichiometry") or []
            net = defaultdict(float)
            for s in stoich:
                net[s["compound"]] += float(s["coefficient"])
            net = {c: v for c, v in net.items() if abs(v) > 1e-12 and c != PROTON}
            parts = [status.get(c, "no structure") for c in net]
            vals = [""] * 4
            if not stoich:
                st = "no stoichiometry"
            elif not net:
                st = "translocation only"
            elif "no structure" in parts:
                st = "missing structure"
            elif "not parsed" in parts or "prediction failed" in parts:
                st = "compound failed"
            elif "R group" in parts:
                st = "R group"
            elif not is_balanced(net, bags):
                st = "unbalanced"
            else:
                vec = np.zeros(n_heads, dtype=np.float64)
                for c, v in net.items():
                    vec += v * heads[c].astype(np.float64)
                m, s = vec.mean(), vec.std()
                if s / KJ_PER_KCAL < DEGENERATE_SIGMA:
                    st = "degenerate"
                else:
                    st = "ok"
                    vals = [f"{m:.3f}", f"{s:.3f}", kcal2(m), kcal2(s)]
            rstatus[st] += 1
            w.writerow([rxn["id"], rxn.get("name", ""), st] + vals)
    os.replace(path + ".tmp", path)
    print("reactions:", dict(rstatus))


if __name__ == "__main__":
    main(_ARGS)
