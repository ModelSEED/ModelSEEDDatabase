```mermaid
%%{init: { 'theme':'dark', 'look':'handDrawn', 'chartOrientation':'horizontal' } }%%
erDiagram

    molecule {
        VARCHAR compound_id PK
        DOUBLE mass
        VARCHAR formula
        TEXT pka
        TEXT pkb
    }

    molecular_structure {
        VARCHAR compound_id FK
        TEXT inchi
        TEXT inchikey
        TEXT smiles
        BINARY molecule
    }

    compartment_annotation {
        VARCHAR reaction_id FK
        VARCHAR ontology_node_id FK
        VARCHAR compartment_index
    }

    catalyst {
        VARCHAR reaction_id PK,FK
        VARCHAR enzyme_id FK
        VARCHAR molecule_id FK
    }
    
    reagent {
        VARCHAR reaction_id PK,FK
        VARCHAR molecule_id PK,FK
        INTEGER compartment_index PK
        DOUBLE stoichiometry
    }

    reaction {
        VARCHAR id
        TEXT name
        BOOLEAN is_transport
    }

    enzyme {
        VARCHAR enzyme_id PK
        VARCHAR enzyme_classification
        VARCHAR description
    }

    reaction ||--o{ reagent : ""
    reagent ||--|| molecule : ""
    reaction ||--o{ catalyst : ""
    catalyst ||--|| enzyme : ""
    catalyst ||--|| molecule : ""
    molecule ||--|| molecular_structure : ""
    reagent ||--|| compartment_annotation : ""

