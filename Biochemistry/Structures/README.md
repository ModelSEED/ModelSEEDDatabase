# ModelSEED Biochemistry Database Structures

Here we keep all the molecular structures we use for the ModelSEED
Biochemistry Database. We download and use molecular structures from
KEGG, MetaCyc and MetaNetX, each stored as SMILES, InChI, and InChIKey
representations in their own folder, and within the
`All_ModelSEED_Structures.txt` file.

We use Marvin to protonate each mol file at a pH of 7 cannot handle
every single file for various reasons, so all files that contain the
word "Original" were generated from the set of uncharged molecular
structure files that were downloaded from their respective sources,
and all files with the word "Charged" concern the files that were
sucessfully protonated by Marvin.

Marvin was also used to generate the pKas and pKbs of the atoms in
each molecular structure.

We also parsed the KEGG molecular structure files for SRUs (structural
repeating units). The information stored for these are passed over by
Marvin when protonating, and not repeated in the resulting output
structure. This caused problems when merging compounds that wre
clearly not the same compound, so we keep a record of which of the
KEGG molecular structure files contain these SRUs.

## Pipeline

The diagram below shows how raw source dumps flow into the
`compound_NN.tsv` files. Each rectangle is a file (or per-source
file family); each rounded box is a script. Read it top-to-bottom:
sources arrive from outside, get derived per-source, get consolidated
across sources, and finally get applied to the canonical compound
TSVs.

```mermaid
%%{init: { 'theme':'dark', 'flowchart': {'curve':'basis'} } }%%
flowchart TD
    subgraph SRC ["Source dumps (one folder per database)"]
        K["KEGG/{InChI,InChIKey,SMILE}_{Original,Charged}Strings.txt"]
        M["MetaCyc/{InChI,InChIKey,SMILE}_{Original,Charged}Strings.txt"]
        C["ChEBI/...Strings.txt"]
        R["Rhea/...Strings.txt"]
        Kp["KEGG/pKa_Strings.txt"]
        Mp["MetaCyc/pKa_Strings.txt"]
        Cp["ChEBI/pKa_Strings.txt (currently unused)"]
        Rp["Rhea/pKa_Strings.txt (currently unused)"]
    end

    subgraph DERIVE ["Per-source derivation"]
        P("Print_Structure_Formula_Charge.py
        — RDKit/OpenBabel parse")
        FC["{KEGG,MetaCyc,ChEBI,Rhea}/*_Formulas_Charges.txt"]
    end
    K & M & C & R --> P --> FC

    subgraph CONS ["Cross-source consolidation"]
        L("List_ModelSEED_Structures.py
        — pick one structure per cpd × type
        priority: Charged-InChI > Original-InChI
        > Charged-SMILE > Original-SMILE
        tie-break: MetaCyc > KEGG > ChEBI > Rhea")
        ALL["All_ModelSEED_Structures.txt"]
        UNI["Unique_ModelSEED_Structures.txt"]
        SR["_reports/Structure_Conflicts.txt
        _reports/Formula_Conflicts.txt
        Pick_Reasons.txt"]
    end
    K & M & C & R & FC --> L
    L --> ALL & UNI & SR

    subgraph OVR ["Curation overrides"]
        ACPS["Curation/overrides/acps_formula_charge.tsv"]
        IGN["Curation/ignores/Ignored_Structures_Publication2020.txt"]
    end

    subgraph APPLY ["Apply to compound_NN.tsv"]
        US("Update_Compound_Structures_Formulas_Charge.py")
        UP("Update_Compound_pKas.py")
        GC("Update_Compound_GroupContribution_Energies.py")
        EQ("Update_Compound_eQuilibrator_Energies.py")
        FCH["compound_NN.tsv: formula, charge, inchikey, smiles"]
        PK["compound_NN.tsv: pka, pkb"]
        DGC["compound_NN.tsv: thermodynamics.GroupContribution"]
        DEQ["compound_NN.tsv: thermodynamics.eQuilibrator
        + overwrites top-level deltag"]
    end
    UNI & ACPS & IGN --> US --> FCH
    Kp & Mp --> UP --> PK
    UNI --> GC
    GTBL["Thermodynamics/ModelSEED/&lcub;KEGG,MetaCyc&rcub;_&lcub;Charged,Original&rcub;_MolAnalysis.tbl"] --> GC --> DGC
    UNI --> EQ
    MNX["Structures/MetaNetX/Structures_in_ModelSEED_and_eQuilibrator.txt"] --> EQ
    ETBL["Thermodynamics/eQuilibrator/MetaNetX_Compound_Energies.tbl"] --> EQ --> DEQ

    REP("Reprint_Biochemistry.py")
    FCH & PK & DGC & DEQ --> REP
    REP --> FINAL["compound_NN.tsv (canonical)
    + compound_NN.json"]

    PROV("Build_Compound_Field_Provenance.py
    — pure data-join, no recompute")
    ALL & UNI & UP & GC & EQ & FINAL --> PROV --> PRV["compound_NN.provenance.tsv
    (per-field source attribution)"]
```

### Where to look when changing things

- **Editing a structure for a single compound.** Find the relevant row
  in the source-DB `*_OriginalStrings.txt` (or the protonated `*_ChargedStrings.txt`)
  for whichever DB owns that alias, then re-run the pipeline. Use
  `Scripts/Structures/Preview_Structure_Update.py <cpd_id> <new_inchi>`
  to see which downstream fields will change before you commit.

- **Understanding why a particular field has its current value.** Look up
  the compound in the matching `compound_NN.provenance.tsv`. Each field
  carries a `<DB>:<external_id>[@stage]` pointer.

- **Following why a particular structure was *picked* in a conflict.** After
  the next pipeline run, `Pick_Reasons.txt` records which branch of the
  picker cascade fired (single_structure, multi_source_agreement,
  smile_db_priority:KEGG, etc).

- **Adding a new protonation engine or pH.** Today this means forking new
  `*_ChargedStrings.txt` files, since the protonation method is encoded
  in the filename. See `sources.yaml` for the schema we are migrating
  toward, where protonation moves into a column.