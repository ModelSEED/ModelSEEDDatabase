use strict;
use fba_tools::fba_toolsImpl;
use Data::Dumper;
use JSON;



sub translation_SSO_to_ModelSEED {

my ($def_hash) = @_;
my $alias_file_url = "../../Biochemistry/Pathways/ModelSEED_Subsystems.tsv";
my $translationfile = "../../Ontologies/KBaseOntology.OntologyTranslation.SSO.ModelSEED.json";

open my $OUTFILE, ">", $translationfile or die "Couldn't open output file $!\n";

open alias_file, $alias_file_url or die "Couldn't open SSO->ModelSEED aliases file $!\n";

my $alias_hash;
my $sso_rxnHash;
my $sso_roleHash;
while (my $inputln = <alias_file>){
    chomp $inputln;
    my @al = split /\t/, $inputln;

    if (!exists $sso_rxnHash->{$al[3]}->{$al[4]}){
        push (@{$alias_hash->{$al[3]}}, $al[4]);
    }
    $sso_rxnHash->{$al[3]}->{$al[4]} = [$al[0],$al[1],$al[2],$al[3],$al[4]];


}


my $EquivalentTerm = {
    equiv_term => "",
    equiv_name => ""
};

my $TransRecord = {
    name => "",
    equiv_terms => []
};

#Replace ontology dictionary references at ontologyX before uploading into workspace
my $OntologyTranslation = {
    comment => "SSO->MS TranlationTable",
    ontology1 => "sso",
    ontology2 => "38284/9/1"
};

my $ontolgoyTable;
foreach my $role (keys ($alias_hash)){

    if (defined $alias_hash->{$role} ){

        my $EquivalentTermArr= [];
        foreach my $r (@{$alias_hash->{$role}}){
            $EquivalentTerm = {
                equiv_term => $r,
                equiv_name => $def_hash->{$r}
            };
            push ($EquivalentTermArr, $EquivalentTerm);
        }

        $TransRecord = {
            name => $role,
            equiv_terms => $EquivalentTermArr
        };
        $OntologyTranslation->{translation}->{$role} = $TransRecord;
    }
}
    my $onTJ = encode_json($OntologyTranslation);
    print $OUTFILE $onTJ;

}


sub translation_KEGG_KO_to_ModelSEED {

my ($def_hash) = @_;
my $alias_file_url = "../../Biochemistry/Aliases/Provenance/KO_modelseed_translations.csv";
my $translationfile = "../../Ontologies/KBaseOntology.OntologyTranslation.KEGG_KO.ModelSEED.json";

open my $OUTFILE, ">", $translationfile or die "Couldn't open output file $!\n";

open alias_file, $alias_file_url or die "Couldn't open KEGG_KO->ModelSEED aliases file $!\n";

    my $alias_hash;
    while (my $inputln = <alias_file>){
        chomp $inputln;
        my @al = split /\t/, $inputln;

        for (my $i =1; $i< @al; $i++){

            push (@{$alias_hash->{$al[0]}}, $al[$i] )
        }
    }

my $Cjson;
    {
        local $/; #Enable 'slurp' mode
        open my $fh, "<", "../../Ontologies/KBaseOntology.OntologyDictionary.KEGG_KO_ontologyDictionary.json";
        $Cjson = <$fh>;
        close $fh;
    }

    my $co = decode_json($Cjson);
    my $cpdStHash;
    my $inchikeyHash;

    my $EquivalentTerm = {
            equiv_term => "",
            equiv_name => ""
        };

    my $TransRecord = {
        name => "",
        equiv_terms => []
    };

    #Replace ontology dictionary references at ontologyX before uploading into workspace
    my $OntologyTranslation = {
        comment => "KEGG_KO->MS TranlationTable",
        ontology1 => "38284/7/7",
        ontology2 => "38284/9/1"
    };


    my $ontolgoyTable;
    foreach my $ko (keys ($co->{term_hash})){

        if (defined $alias_hash->{$ko} ){

            my $EquivalentTermArr= [];
            foreach my $r (@{$alias_hash->{$ko}}){
                $EquivalentTerm = {
                  equiv_term => $r,
                  equiv_name => $def_hash->{$r}
                };
               push ($EquivalentTermArr, $EquivalentTerm);
            }

            $TransRecord = {
              name => $co->{term_hash}->{$ko}->{name},
              equiv_terms => $EquivalentTermArr
            };
            $OntologyTranslation->{translation}->{$ko} = $TransRecord;
        }


    }

    my $onTJ = encode_json($OntologyTranslation);
    print $OUTFILE $onTJ;

}







sub translation_MetaCyc_to_ModelSEED {

my ($def_hash) = @_;
my $alias_file_url = "../../Biochemistry/Aliases/Provenance/Reactions_Aliases.tsv";
my $translationfile = "../../Ontologies/KBaseOntology.OntologyTranslation.MetaCyc_RXN.ModelSEED.json";

open my $OUTFILE, ">", $translationfile or die "Couldn't open output file $!\n";

open alias_file, $alias_file_url or die "Couldn't open Metacyc->ModelSEED aliases file $!\n";

    my $alias_hash;
    while (my $inputln = <alias_file>){
        chomp $inputln;
        my @al = split /\t/, $inputln;
        if ($inputln =~ /KEGG/ ){
            #$alias_hash->{$al[2]} = [$al[0], $al[1],$al[2]];
        }
        elsif ($inputln =~ /MetaCyc/) {
            my @alDot = split /\./, $al[2];
            push (@{$alias_hash->{$alDot[0]}},$al[0] )
            #$alias_hash->{$alDot[0]} = [$al[0], $al[1],$al[2], $alDot[0]];
        }
        else {
            next;
        }
    }

my $Cjson;
    {
        local $/; #Enable 'slurp' mode
        open my $fh, "<", "../../Ontologies/KBaseOntology.OntologyDictionary.MetaCyc_RXN_ontologyDictionary.json";
        $Cjson = <$fh>;
        close $fh;
    }

    my $co = decode_json($Cjson);
    my $cpdStHash;
    my $inchikeyHash;

    my $EquivalentTerm = {
            equiv_term => "",
            equiv_name => ""
        };

    my $TransRecord = {
        name => "",
        equiv_terms => []
    };

    #Replace ontology dictionary references at ontologyX before uploading into workspace
    my $OntologyTranslation = {
        comment => "MetaCyc->MS TranlationTable",
        ontology1 => "38284/10/1",
        ontology2 => "38284/9/1"
    };

    my $ontolgoyTable;
    foreach my $krMETA (keys ($co->{term_hash})){
        my @splitedKr  = split /META:/, $krMETA;
        my $kr= $splitedKr[1];
        #print &Dumper ($alias_hash->{$kr});
        if (defined $alias_hash->{$kr} ){

            # There are duplicate ModelSEED reactions ids(same equation)for the same MetaCyc identifier except for postfix .x
            #rxn02836        RXN-4307.c  MetaCyc
            #rxn23922        RXN-4307.d  MetaCyc
            #rxn23923        RXN-4307.m  MetaCyc

            my $EquivalentTermArr= [];
            foreach my $r (@{$alias_hash->{$kr}}){
                $EquivalentTerm = {
                  equiv_term => $r,
                  equiv_name => $def_hash->{$r}
                };
               push ($EquivalentTermArr, $EquivalentTerm);
            }

            $TransRecord = {
              name => $co->{term_hash}->{$krMETA}->{name},
              equiv_terms => $EquivalentTermArr

            };
            $OntologyTranslation->{translation}->{$kr} = $TransRecord;
        }
    }

    my $onTJ = encode_json($OntologyTranslation);
    print $OUTFILE $onTJ;

}


sub translation_KEGG_rxn_to_ModelSEED {

my ($def_hash) = @_;
my $alias_file_url = "../../Biochemistry/Aliases/Provenance/Reactions_Aliases.tsv";
my $translationfile = "../../Ontologies/KBaseOntology.OntologyTranslation.KEGG_RXN.ModelSEED.json";

open my $OUTFILE, ">", $translationfile or die "Couldn't open output file $!\n";

open alias_file, $alias_file_url or die "Couldn't open MS->KEGG aliases file $!\n";

    my $alias_hash;
    while (my $inputln = <alias_file>){
        chomp $inputln;
        my @al = split /\t/, $inputln;
        if ($inputln =~ /KEGG/ ){
            $alias_hash->{$al[2]} = [$al[0], $al[1],$al[2]];
        }
        elsif ($inputln =~ /MetaCyc/) {
            #$alias_hash->{$al[1]} = [$al[0], $al[1],$al[2]];
        }
        else {
            next;
        }
    }

my $Cjson;
    {
        local $/; #Enable 'slurp' mode
        open my $fh, "<", "../../Ontologies/KBaseOntology.OntologyDictionary.KEGG_RXN_ontologyDictionary.json";
        $Cjson = <$fh>;
        close $fh;
    }

    my $co = decode_json($Cjson);
    my $cpdStHash;
    my $inchikeyHash;

    my $EquivalentTerm = {
	        equiv_term => "",
    	    equiv_name => ""
    	};

    my $TransRecord = {
    	name => "",
    	equiv_terms => []
    };

    #Replace ontology dictionary references at ontologyX before uploading into workspace
    my $OntologyTranslation = {
    	comment => "KEGG->MS TranlationTable",
    	ontology1 => "38284/6/5",
    	ontology2 => "38284/9/1"
    };

    my $ontolgoyTable;
    foreach my $kr (keys ($co->{term_hash})){
    	$EquivalentTerm = [{
	        equiv_term => $alias_hash->{$kr}->[0],
    	    equiv_name => $def_hash->{$kr}->[0]
    	}];

    	$TransRecord = {
    		name => $co->{term_hash}->{$kr}->{name},
    		equiv_terms => $EquivalentTerm

    	};

    	$OntologyTranslation->{translation}->{$kr} = $TransRecord;
    }

    my $onTJ = encode_json($OntologyTranslation);
    print $OUTFILE $onTJ;

}


sub translation_EC_to_ModelSEED {

my ($def_hash) = @_;
my $alias_file_url = "../../Biochemistry/Aliases/Provenance/Enzyme_Class_Reactions_Aliases.tsv";
my $translationfile = "../../Ontologies/KBaseOntology.OntologyTranslation.EBI_EC.ModelSEED.json";
open my $OUTFILE, ">", $translationfile or die "Couldn't open output file $!\n";

open alias_file, $alias_file_url or die "Couldn't open complete data file $!\n";
        my $alias_hash;
        while (my $inputln = <alias_file>){
                chomp $inputln;
                my @al = split /\t/, $inputln;
                $alias_hash->{$al[2]} = $al[0];
        }
    my $Cjson;
    {
        local $/; #Enable 'slurp' mode
        open my $fh, "<", "../../Ontologies/KBaseOntology.OntologyDictionary.EBI_EC_ontologyDictionary.json";
        $Cjson = <$fh>;
        close $fh;
    }
    my $co = decode_json($Cjson);
    my $cpdStHash;
    my $inchikeyHash;

    my $EquivalentTerm = {
          equiv_term => "",
          equiv_name => ""
      };
    my $TransRecord = {
      name => "",
      equiv_terms => []

    };
    #Replace ontology dictionary references at ontologyX before uploading into workspace
    my $OntologyTranslation = {
      comment => "EC->ModelSEED_TranlationTable",
      ontology1 => "38284/8/4",
      ontology2 => "38284/9/1"
    };

    my $ontolgoyTable;
    foreach my $kr (keys ($co->{term_hash})){
        #print "***$alias_hash->{$kr}***\n";
        if (defined $alias_hash->{$kr}){
            my @rxns = split /\|/, $alias_hash->{$kr};


            my $EquivalentTermArr= [];
            foreach my $r (@rxns){
              $EquivalentTerm = {
                  equiv_term => $r,
                  equiv_name => $def_hash->{$r}
              };

              push ($EquivalentTermArr, $EquivalentTerm);
            }
            $TransRecord = {
              name => $co->{term_hash}->{$kr}->{name},
              equiv_terms => $EquivalentTermArr

            };

            $OntologyTranslation->{translation}->{$kr} = $TransRecord;

        }

    }
my $onTJ = encode_json($OntologyTranslation);
print $OUTFILE $onTJ;
}

my $reaction_definitions = "../../Biochemistry/reactions.tsv";
open def_file, $reaction_definitions or die "Couldn't open reaction_definitions file $!\n";

    my $def_hash;
    while (my $inputln = <def_file>){
        chomp $inputln;
        my @al = split /\t/, $inputln;
        $def_hash->{$al[0]} = $al[7];
    }


translation_SSO_to_ModelSEED ($def_hash);
translation_KEGG_KO_to_ModelSEED ($def_hash);
translation_KEGG_rxn_to_ModelSEED ($def_hash);
translation_EC_to_ModelSEED ($def_hash);
translation_MetaCyc_to_ModelSEED ($def_hash);
