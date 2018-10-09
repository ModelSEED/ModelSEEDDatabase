compound_props = {
    "id": {"type": "string", "pattern": "^cpd\d{5}$"},
    "abbreviation": {"type": "string"},
    "name": {"type": "string"},
    "formula": {"type": "string"},
    "source": {"type": "string"},
    "inchikey": {"type": "string"},
    "smiles": {"type": "string"},
    "aliases": {"type": "string"},
    "linked_compound": {"type": "string", "pattern": "^cpd\d{5}$|^null$"},
    "comprised_of": {"type": "string"},
    "charge": {"type": "string"},
    "mass": {"type": "string"},
    "deltag": {"type": ["number", "string"]},
    "deltagerr": {"type": ["number", "string"]},
    "pka": {"type": ["number", "string"]},
    "pkb": {"type": ["number", "string"]},
    "is_core": {"type": "number"},
    "is_obsolete": {"type": "number"},
    "is_cofactor": {"type": "number"},
    "abstract_compound": {"type": "string"},
}

compounds = {
    "type": "object",
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
    "aliases": {"type": "string"},
    "linked_reaction": {"type": "string", "pattern": "(rxn\d{5})*"},
    "direction": {"type": "string", "pattern": "^[<=>]$"},
    "reversibility": {"type": "string", "pattern": "^[<=>?]$"},
    "ec_numbers": {"type": "string"},
    "equation": {"type": "string"},
    "deltag": {"type": ["number", "string"]},
    "deltagerr": {"type": ["number", "string"]},
    "notes": {"type": "string"},
    "pathways": {"type": "string"},
    "status": {"type": "string"},
    "stoichiometry": {"type": "string"},
    "is_transport": {"type": "number"},
    "is_obsolete": {"type": "number"},
    "abstract_reaction": {"type": "string"},
    "source": {"type": "string"},
}

reactions = {
    "type": "object",
    "patternProperties": {
        "rxn\d{5}": {
            "type": "object",
            "properties":  reaction_props,
            "required": list(reaction_props.keys()),
            "additionalProperties": False,
        }
    }
}
