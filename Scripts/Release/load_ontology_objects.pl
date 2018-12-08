use strict;
use fba_tools::fba_toolsImpl;
use Data::Dumper;
use JSON;

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
my $wsClient = Bio::KBase::kbaseenv::ws_client();
#my $workspace_name = "KBaseOntology";
my $workspace_name = "janakakbase:narrative_1543604022660"; #test prod ws
#my $workspace_name = "janakakbase:1466175040023"; #test appdev ws

system ("python create_ontology_dictionaries.py");
system ("perl build_ontology_translation_objects.pl");

opendir my $dir, "../../Ontologies/" or die "Cannot open Ontologies directory: $!";
my @files = readdir $dir;
closedir $dir;
print &Dumper (\@files);
foreach my $f (@files){

	if ($f =~ /KBaseOntology.OntologyDictionary/){

		my @splitN = split /\./, $f;
		my $Ojson;
    	{
        local $/; #Enable 'slurp' mode
        open my $fh, "<", "../../Ontologies/$f";
        $Ojson = <$fh>;
        chomp $Ojson;
        close $fh;
    	}
        my $co = decode_json($Ojson);

		my $obj_info_list = undef;
    	eval {
	        $obj_info_list = $wsClient->save_objects({
	            'workspace'=>$workspace_name,
	            'objects'=>[{
	            'type'=> $splitN[0].".".$splitN[1],
	            'data'=>$co,
	            'name'=>$splitN[2]
	            }]
	        });
    	};
    	if ($@) {
        	die "Error saving modified Ontology object to workspace:\n".$@;
    	}

        print &Dumper ($obj_info_list);
	}


}
