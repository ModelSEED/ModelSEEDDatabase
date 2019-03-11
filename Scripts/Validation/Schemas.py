compound_props = {
    "id": {"type": "string", "pattern": "^cpd\d{5}$"},
    "abbreviation": {"type": "string"},
    "name": {"type": "string"},
    "formula": {"type": ["string","null"]},
    "source": {"type": "string"},
    "inchikey": {"type": "string"},
    "smiles": {"type": "string"},
    "aliases": {"type": ["string","null"]},
    "linked_compound": {"type": ["string","null"], "pattern": "^cpd\d{5}$|^null$"},
    "comprised_of": {"type": ["string","null"]},
    "charge": {"type": "string"},
    "mass": {"type": ["string","null"]},
    "deltag": {"type": ["number","string","null"]},
    "deltagerr": {"type": ["number","string","null"]},
    "pka": {"type": ["number", "string"]},
    "pkb": {"type": ["number", "string"]},
    "is_core": {"type": "number"},
    "is_obsolete": {"type": "number"},
    "is_cofactor": {"type": "number"},
    "abstract_compound": {"type": ["string","null"]},
}

compounds = {
    "type": "array",
    "patternProperties": {
        "cpd\d{5}": {
            "type": "object",
            "properties":  compound_props,
            "required": list(compound_props.keys()),
            "additionalProperties": False,
        }
    },
    "additionalProperties": False,
}

reaction_props = {
    "id": {"type": "string", "pattern": "^rxn\d{5}$"},
    "abbreviation": {"type": "string"},
    "name": {"type": "string"},
    "code": {"type": "string"},
    "compound_ids": {"type": "string"},
    "definition": {"type": "string"},
    "aliases": {"type": ["string","null"]},
    "linked_reaction": {"type": ["string","null"], "pattern": "(rxn\d{5})*"},
    "direction": {"type": "string", "pattern": "^[<=>]$"},
    "reversibility": {"type": "string", "pattern": "^[<=>?]$"},
    "ec_numbers": {"type": ["string","null"]},
    "equation": {"type": "string"},
    "deltag": {"type": ["number","string","null"]},
    "deltagerr": {"type": ["number","string","null"]},
    "notes": {"type": ["string","null"]},
    "pathways": {"type": ["string","null"]},
    "status": {"type": "string"},
    "stoichiometry": {"type": "string"},
    "is_transport": {"type": "number"},
    "is_obsolete": {"type": "number"},
    "abstract_reaction": {"type": ["string","null"]},
    "source": {"type": "string"},
}

reactions = {
    "type": "array",
    "patternProperties": {
        "rxn\d{5}": {
            "type": "object",
            "properties":  reaction_props,
            "required": list(reaction_props.keys()),
            "additionalProperties": False,
        }
    }
}
