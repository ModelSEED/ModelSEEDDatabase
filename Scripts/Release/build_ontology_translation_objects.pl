use strict;
use fba_tools::fba_toolsImpl;
use Data::Dumper;
use JSON;


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
    	ontology1 => "KEGG",
    	ontology2 => "ModelSEED"
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
      ontology1 => "EC",
      ontology2 => "ModelSEED"
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

translation_KEGG_rxn_to_ModelSEED ($def_hash);
translation_EC_to_ModelSEED ($def_hash);
