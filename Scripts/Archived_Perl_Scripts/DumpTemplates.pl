#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );
$|=1;

(my $data,my $info) = get_workspace_object("KBaseTemplateModels/GramNegModelTemplate");

my $directory = $Bin."/../Templates/GramNegative/";
print $directory."\n";

my $biocpds = $data->{templateBiomasses};
open ( BIOMASS, ">", $directory."/BiomassCompounds.txt");
print BIOMASS "cpdid\tcompartment\tclass\tcoefficient\tcoefficient type\tlinked compounds\n";
for (my $i=0; $i < @{$biocpds}; $i++) {
	my $line = "";
	if ($biocpds->[$i]->{compound_ref} =~ m/\/([^\/]+)$/) {
		$line .= $1;
	}
	$line .= "\t";
	if ($biocpds->[$i]->{compartment_ref} =~ m/\/([^\/]+)$/) {
		$line .= $1;
	}
	$line .= "\t".$biocpds->[$i]->{class}."\t".$biocpds->[$i]->{coefficient}."\t".$biocpds->[$i]->{coefficientType}."\t";
	if (defined($biocpds->[$i]->{linked_compound_refs})) {
		my $links = $biocpds->[$i]->{linked_compound_refs};
		for (my $j=0; $j < @{$links}; $j++) {
			if ($links->[$j] =~ m/\/([^\/]+)$/) {
				if ($j > 0) {
					$line .= "|";
				}
				$line .= $1.":".$biocpds->[$i]->{link_coefficients}->[$j];
			}
		}
	}
	print BIOMASS $line."\n";
}
close(BIOMASS);

my $rxns = $data->{templateReactions};
open ( TEMPLATE, ">", $directory."/Reactions.txt");
print TEMPLATE "rxnid\tcompartment\tdirection\tgfdir\ttype\tbasecost\tforpen\trevpen\tcomplexes\n";
for (my $i=0; $i < @{$rxns}; $i++) {
	my $line = "";
	if ($rxns->[$i]->{reaction_ref} =~ m/\/([^\/]+)$/) {
		$line .= $1;
	}
	$line .= "\t";
	if ($rxns->[$i]->{compartment_ref} =~ m/\/([^\/]+)$/) {
		$line .= $1;
	}
	if (!defined()) {
		$rxns->[$i]->{GapfillDirection} = "=";
	}
	$line .= "\t".$rxns->[$i]->{direction}."\t".$rxns->[$i]->{GapfillDirection}."\t".$rxns->[$i]->{type}."\t".$rxns->[$i]->{base_cost}."\t".$rxns->[$i]->{forward_penalty}."\t".$rxns->[$i]->{reverse_penalty}."\t";
	if (defined($rxns->[$i]->{complex_refs})) {
	my $links = $rxns->[$i]->{complex_refs};
	for (my $j=0; $j < @{$links}; $j++) {
		if ($links->[$j] =~ m/\/([^\/]+)$/) {
			if ($j > 0) {
				$line .= "|";
			}
			$line .= $1;
		}
	} 
	}
	print TEMPLATE $line."\n";
}
close(TEMPLATE);