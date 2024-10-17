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

    molecule_annotation {
        VARCHAR molecule_id FK
        VARCHAR ontology_node_id FK
    }

    molecular_structure {
        VARCHAR molecule_id FK
        TEXT inchi
        TEXT inchikey
        TEXT smiles
        BINARY molecule
    }

    enzyme_abstraction {
        VARCHAR enzyme_id FK
        VARCHAR ontology_node_id FK
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

    reaction_annotation {
        VARCHAR reaction_id FK
        VARCHAR ontology_node_id FK
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
    molecule ||--o{ molecule_annotation : ""
    reaction ||--o{ reaction_annotation : ""
    enzyme ||--o{ enzyme_abstraction : ""
    reagent ||--|| compartment_annotation : ""

