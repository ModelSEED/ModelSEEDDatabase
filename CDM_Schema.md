```mermaid
%%{init: { 'theme':'dark', 'look':'handDrawn', 'chartOrientation':'horizontal' } }%%
erDiagram
    compound {
        VARCHAR compound_id PK
        VARCHAR structure_id FK
        DOUBLE mass
        VARCHAR formula
        DOUBLE charge
        TEXT pka
        TEXT pkb
    }

    molecular_structure {
        VARCHAR id
        VARCHAR compound_id FK
        TEXT inchi
        TEXT inchikey
        TEXT smiles
        BINARY molecule
    }

    entity_as_alias {
        VARCHAR alias_id FK
        VARCHAR compound_id FK
        VARCHAR reaction_id FK
        VARCHAR pathway_id FK
        VARCHAR molecular_structure_id FK
    }

    alias {
        VARCHAR id
        VARCHAR database
        VARCHAR identifier
    }
    
    reagent {
        VARCHAR reaction_id PK,FK
        VARCHAR compound_id PK,FK
        INTEGER compartment_index PK
        DOUBLE stoichiometry
    }

    reaction {
        VARCHAR id
        TEXT name
        BOOLEAN is_transport
    }

    reaction_in_pathway {
        VARCHAR reaction_id FK
        VARCHAR pathway_id FK
    }

    pathway {
        VARCHAR id
        TEXT name
        TEXT description
    }

    pathway ||--o{ reaction_in_pathway : "contains"
    reaction ||--o{ reaction_in_pathway : "contains"
    reaction ||--o{ reagent : "catalyzes"
    reagent ||--|| compound : "involves"
    compound ||--|| molecular_structure : "defined by"
    compound ||--o{ entity_as_alias : "represented as"
    reaction ||--o{ entity_as_alias : "represented as"
    pathway ||--o{ entity_as_alias : "represented as"
    molecular_structure ||--o{ entity_as_alias : "represented as"
    alias ||--o{ entity_as_alias : "represented as"
