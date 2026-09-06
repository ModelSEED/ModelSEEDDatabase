# ModelSEED Biochemistry Database — Thermodynamics data

Per-source thermodynamic estimates for every compound and reaction in
the database, stored additively so no source is ever silently overwritten.

The estimates come from four sources today:

- **Group Contribution** — MFAToolkit's implementation of the
  Jankowski et al. 2008 method (<a href="https://doi.org/10.1529/biophysj.107.124784">10.1529/biophysj.107.124784</a>),
  computed against each compound's current SMILES.
- **eQuilibrator** — Noor et al. 2013 (<a href="https://doi.org/10.1371/journal.pcbi.1003098">10.1371/journal.pcbi.1003098</a>);
  values pulled from eQuilibrator's public cache via InChIKey matching.
- **dGPredictor** — Wang et al. 2021 (<a href="https://doi.org/10.1371/journal.pcbi.1009448">10.1371/journal.pcbi.1009448</a>);
  a group-decomposition + ML model.
- **dGPredictor-ModelSEED** — the dGPredictor model retrained on the
  current ModelSEED compound corpus.

Values are kept as per-source dicts under the `thermodynamics` key on
every compound and reaction JSON record. The dict shape is:

```json
"thermodynamics": {
    "Group contribution":    [4.15, 1.22, "="],
    "eQuilibrator":          [-3.46, 0.05, ">"],
    "dGPredictor":           [-3.82, 0.02, ">"],
    "dGPredictor-ModelSEED": [-3.77, 0.87, ">"]
}
```

Reactions carry `[energy, error, operator]` triples where the operator
(`>`, `<`, `=`, or `?`) is that source's own thermodynamic-direction call.
Compounds carry `[energy, error]` pairs (no direction — compounds don't
have a reaction direction).

The sentinel value `10000000` in the energy or error slot means "no
usable estimate from this source"; operator `?` means "unknown direction."

---

## Convention A: pH 7 in Chris Henry's 2010 formalism

**Every `thermodynamics["Group contribution"]` entry ships in
Convention A.** This is a load-bearing property of the database — mixing
conventions within a reaction gives wrong Δ<sub>r</sub>G' by
(Δn<sub>H</sub> − n<sub>transported H⁺</sub>) × 9.539 kcal/mol.

**Convention A** (as authored in the 2010 ModelSEED biochemistry paper
and produced by MFAToolkit natively):

1. Compound Δ<sub>f</sub>G *includes* full hydrogen accounting. Water
   is −56.687 kcal/mol (accounting for its 2 H atoms in elemental
   formation), not −37.6 kcal/mol.
2. Free proton is assigned a fixed formation energy of
   **Δ<sub>f</sub>G(H⁺) = −9.5 kcal/mol**. This value equals
   −RT ln(10) × 7 at 298 K — the pH 7 correction lives here.
3. No per-compound Legendre transform. When a reaction is balanced, any
   H⁺ terms carry the pH 7 correction automatically through the H⁺
   formation energy.

### Convention A vs Convention B — same physics, different bookkeeping

Alberty's convention (which eQuilibrator emits) is mathematically
equivalent to Convention A at pH 7 but keeps the books differently:
per-compound Δ<sub>f</sub>G' has the H contribution subtracted at the
compound level (water = −37.6 kcal/mol), and H⁺ has Δ<sub>f</sub>G' = 0.

Both conventions produce the same Δ<sub>r</sub>G' for any properly-
balanced reaction. But you cannot mix them within a single reaction.

### Worked example: water dissociation (H₂O → OH⁻ + H⁺)

**Convention A:**
- Δ<sub>f</sub>G(H₂O) = −56.687
- Δ<sub>f</sub>G(OH⁻) = ~ −37.6
- Δ<sub>f</sub>G(H⁺) = −9.5
- Δ<sub>r</sub>G' = −37.6 + (−9.5) − (−56.687) = **+9.587 kcal/mol**

**Convention B (Alberty):**
- Δ<sub>f</sub>G'(H₂O) = −37.6
- Δ<sub>f</sub>G'(OH⁻) = ~ −18.9  (i.e. −37.6 − 1×(−9.539) with 1H shift)
- Δ<sub>f</sub>G'(H⁺) = 0
- Δ<sub>r</sub>G' = −18.9 + 0 − (−37.6) = **+9.7 kcal/mol**

Same physics. Different accounting.

### Load-bearing Convention A anchors

Nine compounds carry Convention A values injected at the compound level
by `Scripts/Thermodynamics/Update_Compound_GroupContribution_Energies.py`
rather than derived from MFAToolkit's per-source MolAnalysis alias
lookup. MFAToolkit's atom-labeler can't decompose these standalone (its
cues only match them as substructures of larger molecules), so the
alias-lookup returns sentinel. The injected values come from Chris
Henry's own cue-database ENERGY entries and have been independently
cross-validated against eQuilibrator (Convention B) values shifted to
Convention A via `Δ<sub>f</sub>G<sub>A</sub> = Δ<sub>f</sub>G<sub>B</sub> − n<sub>H</sub> × 9.539`.
Eight agree to within 0.5 kcal/mol; the ninth (H₂S) has a 2.88 kcal/mol
residual explained by eQuilibrator's pseudo-species treatment vs
Chris's neutral-form value.

| Compound | Formula | Δ<sub>f</sub>G (kcal/mol) | Notes |
|---|---|---:|---|
| cpd00067 (H⁺) | H | −9.5 | Chris's pH 7 correction anchor |
| cpd00001 (H₂O) | H₂O | −56.687 | Ubiquitous; unlocks GC coverage on ~20K reactions |
| cpd00013 (NH₄⁺) | H₄N | −18.97 | |
| cpd00011 (CO₂) | CO₂ | −92.26 | |
| cpd00242 (HCO₃⁻) | CHO₃ | −140.26 | |
| cpd00025 (H₂O₂) | H₂O₂ | −32.05 | |
| cpd11640 (H₂) | H₂ | 4.2065 | |
| cpd00007 (O₂) | O₂ | 3.9197 | |
| cpd00239 (H₂S) | HS | −6.66 | Formula shows the Charged (pH 7) HS⁻ form; identity is H₂S. eQuilibrator's pseudo-species value gives a 2.88 kcal/mol offset via the same math. |

The `Group contribution` entry on every other compound comes from the
alias-lookup pipeline over Chris's per-source mol-file corpora (see the
`ModelSEED/*_MolAnalysis.tbl` files in this directory and
`Scripts/Thermodynamics/Update_Compound_GroupContribution_Energies.py`
for the mean-across-aliases resolver).

---

## Non-GC sources still ship in Convention B

`thermodynamics["eQuilibrator"]`, `thermodynamics["dGPredictor"]`, and
`thermodynamics["dGPredictor-ModelSEED"]` currently ship in **Convention B**
(what those tools emit). Bringing them into Convention A is a separate
follow-up PR that requires per-compound n<sub>H</sub> lookup and a
uniform transform pass.

Downstream tools consuming these sources should be aware of this mixture.
A convention-aware wrapper: if `label == "Group contribution"` use as-is
with H⁺ = −9.5 kcal/mol; for the other labels use as-is with H⁺ = 0
kcal/mol.

---

## Top-level `deltag` / `deltagerr` / `reversibility` — historical scalars

Each compound and reaction also carries top-level scalar fields
(`deltag`, `deltagerr`, and — for reactions — `reversibility`) picked
from one of the per-source entries by
`Scripts/Thermodynamics/Promote_Reaction_Thermodynamics_to_Canonical.py`
using a tier + lowest-uncertainty policy.

**These top-level scalars are historical values with mixed conventions**
(a mix of Convention A and Convention B, depending on which source was
picked for each record) and are **preserved as-is by any GC-cleanup
work** for backwards compatibility.

**⚠ DEPRECATION NOTICE ⚠**

Top-level `deltag`, `deltagerr`, and `reversibility` are scheduled for
retirement in a future release once consumers have migrated to reading
the per-source dict. The migration path is:

```python
# Old (deprecated — historical mixed-convention scalar):
dg = reaction['deltag']

# New (preferred — per-source, Convention A for Group contribution):
dg = reaction['thermodynamics']['Group contribution'][0]
```

If your tool consumes `deltag` directly today, please move to the
`thermodynamics[<source>]` accessor at your convenience — the top-level
fields will be removed in a future release. The exact release cadence
will be announced separately.

---

## Reproducibility: pinned MFAToolkit version

The `Group contribution` values are produced by the MFAToolkit tool
(<a href="https://github.com/ModelSEED/MFAToolkit">github.com/ModelSEED/MFAToolkit</a>),
run over the current ModelSEED compound SMILES. The specific tool
version used to produce the current values is pinned in

    ModelSEED/MFAToolkit_version.txt

and the raw tool output — a compound-level `MolAnalysis.tbl` — is
committed in this repository under

    ModelSEED/<name>_MolAnalysis.tbl

so anyone with the ModelSEEDDatabase clone can (a) see which MFAToolkit
release produced the current values, (b) install that release, (c)
rerun it against the current compound SMILES, and (d) diff the result.

To regenerate the `Group contribution` values yourself:

```bash
# 1. Install the pinned MFAToolkit release
git clone https://github.com/ModelSEED/MFAToolkit
cd MFAToolkit && git checkout $(cat ../ModelSEEDDatabase/Biochemistry/Thermodynamics/ModelSEED/MFAToolkit_version.txt)
make mfatoolkit

# 2. Run over the current ModelSEED SMILES
# (script in Scripts/Thermodynamics/ produces the SDF, invokes mfatoolkit,
#  and updates thermodynamics["Group contribution"] on every compound)
cd ../ModelSEEDDatabase/Scripts/Thermodynamics
./Rerun_GroupContribution.sh
```

---

## See also

- `Scripts/Thermodynamics/README.md` — the update-pipeline scripts,
  including which script writes which key.
- `Papers/NAR_Update_2026/data/gc_drift_source_2026-08-06.md` — the
  investigation that identified the 42/58 Convention-A/B mix in the
  historical top-level `deltag` and the fix path taken here.
- `Papers/NAR_Update_2026/data/gc_cleanup_pr_plan_2026-08-06.md` — the
  cleanup PR plan this README ships as part of.
- Jankowski et al. 2008 (<a href="https://doi.org/10.1529/biophysj.107.124784">10.1529/biophysj.107.124784</a>) — the original Group Contribution method for biochemistry.
- Alberty, R.A. *Thermodynamics of Biochemical Reactions* (2003) — the
  canonical reference for the transformed (Convention B) formalism.
