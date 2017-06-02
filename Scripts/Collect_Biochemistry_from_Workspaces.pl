#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header=1;

#######################################################
#Initialization
#######################################################

use Bio::KBase::fbaModelServices::ScriptHelpers qw( getToken );
my $AToken = getToken();

use Bio::KBase::fbaModelServices::Impl;
my $FBAImpl = Bio::KBase::fbaModelServices::Impl->new({'fbajobcache' => "/homes/seaver/Projects/KBase_Scripts/FBA_Scripts/JobsCache",
                                                       'jobserver-url' => "http://kbase.us/services/workspace",
                                                       'fbajobdir' => "/tmp/fbajobs",
                                                       'mfatoolkitbin' => "/homes/chenry/kbase/MFAToolkit/bin/mfatoolkit",
                                                       'probanno-url' => "http://140.221.85.86:7073/",
                                                       'mssserver-url' => "http://bio-data-1.mcs.anl.gov/services/ms_fba",
                                                       'accounttype' => "kbase",
                                                       'workspace-url' => "http://kbase.us/services/ws",
                                                       'defaultJobState' => "queued",
                                                       'gaserver-url' => "http://kbase.us/services/genome_annotation",
                                                       'idserver-url' => "http://kbase.us/services/idserver"});
$FBAImpl->_setContext(undef,{auth=>$AToken});
my $WSClient = $FBAImpl->_workspaceServices();

my %Rxns=();
my %Cpds=();

foreach my $template (@{$WSClient->list_objects({workspaces=>["KBaseTemplateModels"],type=>"KBaseFBA.ModelTemplate"})}){
    my $templateObj = $FBAImpl->_get_msobject("ModelTemplate","KBaseTemplateModels",$template->[1]);

    foreach my $tmplrxn (@{$templateObj->templateReactions()}){
	next unless $tmplrxn->type() eq "universal" || $tmplrxn->type() eq "conditional";

	my $rxn = $tmplrxn->reaction();
	$Rxns{$rxn->id()}{$template->[1]}=1;

	foreach my $rgt (@{$rxn->reagents()}){
	    $Cpds{$rgt->compound()->id()}{$template->[1]}=1;
	}
    }
}

open(OUT, "> KBaseTemplateModels.cpd");
print OUT join("\n",map { $_."\t".join("|",sort keys %{$Cpds{$_}}) } sort keys %Cpds),"\n";
close(OUT);

open(OUT, "> KBaseTemplateModels.rxn");
print OUT join("\n",map { $_."\t".join("|",sort keys %{$Rxns{$_}}) } sort keys %Rxns),"\n";
close(OUT);

__END__
undef(%Cpds);
undef(%Rxns);

foreach my $media (@{$WSClient->list_objects({workspaces=>["KBaseMedia"],type=>"KBaseBiochem.Media"})}){
    my $mediaObj = $FBAImpl->_get_msobject("Media","KBaseMedia",$media->[1]);

    foreach my $medcpd (@{$mediaObj->mediacompounds()}){
	$Cpds{$medcpd->compound()->id()}=1;
    }
}

open(OUT, "> KBaseMedia.cpd");
print OUT join("\n", sort keys %Cpds),"\n";
close(OUT);

foreach my $ws ("KBasePublicModelsV4","PlantSEED"){
    foreach my $model (@{$WSClient->list_objects({workspaces=>[$ws],type=>"KBaseFBA.FBAModel"})}){
	undef(%Cpds);
	undef(%Rxns);
	
	my $modelObj = $FBAImpl->_get_msobject("FBAModel",$ws,$model->[1]);
       
	foreach my $mdlcpd (@{$modelObj->modelcompounds()}){
	    $Cpds{$mdlcpd->compound()->id()}=1;
	}
	foreach my $mdlrxn (@{$modelObj->modelreactions()}){
	    $Rxns{$mdlrxn->reaction()->id()}=1;
	}
    }
    
    open(OUT, "> ".$ws.".cpd");
    print OUT join("\n", sort keys %Cpds),"\n";
    close(OUT);
	
    open(OUT, "> ".$ws.".rxn");
    print OUT join("\n", sort keys %Rxns),"\n";
    close(OUT);
}
