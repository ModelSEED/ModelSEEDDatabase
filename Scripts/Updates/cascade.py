"""
Pipeline cascade: maps each manifest type to the ordered list of
pipeline stages that need to re-run after the manifest is applied.

Each stage is a (label, command) tuple. Apply_Manifest prints the
labels; Refresh_Pipeline.py executes the commands.
"""
import os

SCRIPTS = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))


def _cmd(rel_path):
    """Build a runnable command from a path relative to Scripts/."""
    return ['python3', os.path.join(SCRIPTS, rel_path)]


# Each stage: stable label + command to run.
STAGES = {
    'Print': (
        'Print_Structure_Formula_Charge',
        _cmd('Structures/Print_Structure_Formula_Charge.py'),
    ),
    'List': (
        'List_ModelSEED_Structures',
        _cmd('Structures/List_ModelSEED_Structures.py'),
    ),
    'UpdateStructures': (
        'Update_Compound_Structures_Formulas_Charge',
        _cmd('Structures/Update_Compound_Structures_Formulas_Charge.py'),
    ),
    'UpdatePkas': (
        'Update_Compound_pKas',
        _cmd('Structures/Update_Compound_pKas.py'),
    ),
    'AddPkaSources': (
        'Add_Compound_pKa_Sources',
        _cmd('Structures/Add_Compound_pKa_Sources.py'),
    ),
    'UpdateGC': (
        'Update_Compound_GroupContribution_Energies',
        _cmd('Thermodynamics/Update_Compound_GroupContribution_Energies.py'),
    ),
    'UpdateEQ': (
        'Update_Compound_eQuilibrator_Energies',
        _cmd('Thermodynamics/Update_Compound_eQuilibrator_Energies.py'),
    ),
    'UpdateAliases': (
        'Update_Compound_Aliases',
        _cmd('Biochemistry/Refresh_DB_after_Changes/Update_Compound_Aliases.py'),
    ),
    'Reprint': (
        'Reprint_Biochemistry',
        _cmd('Biochemistry/Reprint_Biochemistry.py'),
    ),
    'BuildProvenance': (
        'Build_Compound_Field_Provenance',
        _cmd('Provenance/Build_Compound_Field_Provenance.py'),
    ),
    'ValidateFAISS': (
        'Validate_FAISS_Outputs (Tier 1)',
        _cmd('Validation/Validate_FAISS_Outputs.py') + ['--mode=smiles'],
    ),
}


# Per-manifest-type cascade in execution order. The runner skips
# unknown stages with a warning.
CASCADE = {
    'structure_update': [
        'Print', 'List', 'UpdateStructures', 'UpdatePkas', 'AddPkaSources',
        'UpdateGC', 'UpdateEQ', 'Reprint', 'BuildProvenance',
        'ValidateFAISS',
    ],
    'protonation_replace': [
        'List', 'UpdateStructures', 'UpdatePkas', 'AddPkaSources',
        'UpdateGC', 'UpdateEQ', 'Reprint', 'BuildProvenance',
        'ValidateFAISS',
    ],
    'override_add': [
        'UpdateStructures', 'Reprint', 'BuildProvenance',
    ],
    'override_remove': [
        'UpdateStructures', 'Reprint', 'BuildProvenance',
    ],
    'ignore_add': [
        'List', 'UpdateStructures', 'Reprint', 'BuildProvenance',
    ],
    'alias_add': [
        'UpdateAliases', 'BuildProvenance',
    ],
    'alias_remove': [
        'UpdateAliases', 'BuildProvenance',
    ],
    'pka_replace': [
        'UpdatePkas', 'AddPkaSources', 'Reprint', 'BuildProvenance',
    ],
}


def stages_for(manifest_type):
    """Return the ordered list of stage keys for a manifest type."""
    if manifest_type not in CASCADE:
        raise KeyError(f"unknown manifest type: {manifest_type!r}; "
                       f"known: {sorted(CASCADE.keys())}")
    return CASCADE[manifest_type]
