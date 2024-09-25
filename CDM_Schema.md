```mermaid
%%{init: { 'theme':'dark', 'look':'handDrawn', 'chartOrientation':'horizontal' } }%%
erDiagram
    compound {
        VARCHAR compound_id PK
        ARRAYTYPE structures
        ARRAYTYPE aliases
        ARRAYTYPE reagents
        DOUBLE mass
        VARCHAR formula
        DOUBLE charge
        TEXT pka
        TEXT pkb
    }

    molecular_structure {
        VARCHAR structure_id
        ARRAYTYPE compound_ids FK
        TEXT inchi
        TEXT inchikey
        TEXT smiles
        BINARY molecule
    }

    alias {
        VARCHAR alias_id
        ARRAYTYPE compound_ids
        VARCHAR database
        VARCHAR identifier
    }
    
    reagent {
        VARCHAR reagent_id
        TEXT compound_id
        INTEGER compartment_index
        DOUBLE stoichiometry
    }

    reaction {
        VARCHAR reaction_id
        TEXT name
        ARRAYTYPE reagents
        ARRAYTYPE aliases
        ARRAYTYPE pathways
        BOOLEAN is_transport
    }

    pathway {
        VARCHAR pathway_id
        TEXT name
        ARRAYTYPE aliases
        ARRAYTYPE reactions
    }

    pathway ||--|{ reaction : "contains"
    reaction ||--o{ reagent : "catalyzes"
    reagent ||--|| compound : "involves"
    compound ||--|{ molecular_structure : "defined by"
    compound |{--|{ alias : "represented as"
