#!/usr/bin/env perl
use warnings;
use strict;

#######################################################
#Initialization
#######################################################

use Bio::KBase::fbaModelServices::ScriptHelpers qw( getToken );
my $AToken = getToken();

use Bio::KBase::fbaModelServices::Impl;
my $FBAImpl = Bio::KBase::fbaModelServices::Impl->new({'fbajobcache' => "/homes/seaver/Projects/KBase_Scripts/FBA_Scripts/JobsCache",
						       'jobserver-url' => "http://kbase.us/services/workspace",
						       'fbajobdir' => "/tmp/fbajobs",
						       'mfatoolkitbin' => "/vol/model-prod/kbase/MFAToolkit/bin/mfatoolkit",
#						       'mfatoolkitbin' => "/homes/seaver/Software/MFAToolkit/bin/mfatoolkit",
						       'probanno-url' => "http://140.221.85.86:7073/",
						       'mssserver-url' => "http://bio-data-1.mcs.anl.gov/services/ms_fba",
						       'accounttype' => "kbase",
						       'workspace-url' => "http://kbase.us/services/ws",
						       'defaultJobState' => "queued",
						       'gaserver-url' => "http://kbase.us/services/genome_annotation",
						       'idserver-url' => "http://kbase.us/services/idserver"});
$FBAImpl->_setContext(undef,{auth=>$AToken});

my $bioObj = $FBAImpl->_get_msobject("Biochemistry","kbase","plantdefault");

my %Rxns_Pathways=();
my %Pathways=();
foreach my $rxnset (@{$bioObj->reactionSets()}){
    $Pathways{$rxnset->type()}=1;

    foreach my $ref (@{$rxnset->reaction_refs()}){
	my $obj = $bioObj->getLinkedObject($ref);
	$Rxns_Pathways{$obj->id()}{$rxnset->type()}{$rxnset->id()}=1;
    }
}

open(OUT, "> ../Pathways/plantdefault.pathways.tsv");
print OUT "ID\t".join("\t",sort keys %Pathways)."\n";
foreach my $rxn (sort { $a cmp $b } keys %Rxns_Pathways){
    print OUT $rxn,"\t",join("\t", map { exists($Rxns_Pathways{$rxn}{$_}) ? join("|", sort keys %{$Rxns_Pathways{$rxn}{$_}}) : "null" } sort keys %Pathways ),"\n"; 
}
close(OUT);
