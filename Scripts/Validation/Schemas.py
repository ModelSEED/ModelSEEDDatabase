compound_props = {
    "id": {"type": "string", "pattern": "^cpd\d{5}$"},
    "abbreviation": {"type": "string"},
    "name": {"type": "string"},
    "formula": {"type": "string"},
    "source": {"type": "string"},
    "structure": {"type": "string"},
    "aliases": {"type": "string"},
    "linked_compound": {"type": "string", "pattern": "^cpd\d{5}$|^null$"},
    "comprised_of": {"type": "string"},
    "charge": {"type": "string"},
    "mass": {"type": "string"},
    "deltag": {"type": "string"},
    "deltagerr": {"type": "string"},
    "pka": {"type": "string"},
    "pkb": {"type": "string"},
    "is_core": {"type": "string"},
    "is_obsolete": {"type": "string"},
    "is_cofactor": {"type": "string"},
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
    "id": {"type": "string", "pattern": "cpd\d{5}"},
    "abbreviation": {"type": "string"},
    "name": {"type": "string"},
    "formula": {"type": "string"},
    "source": {"type": "string"},
    "structure": {"type": "string"},
    "aliases": {"type": "string"},
    "linked_compound": {"type": "string"},
    "comprised_of": {"type": "string"},
    "charge": {"type": "string"},
    "mass": {"type": "string"},
    "deltag": {"type": "string"},
    "deltagerr": {"type": "string"},
    "pka": {"type": "string"},
    "pkb": {"type": "string"},
    "is_core": {"type": "string"},
    "is_obsolete": {"type": "string"},
    "is_cofactor": {"type": "string"},
    "abstract_compound": {"type": "string"},
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
