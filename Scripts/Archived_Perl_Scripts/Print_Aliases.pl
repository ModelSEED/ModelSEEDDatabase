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


#######################################################
#Biochemistry and Model
#######################################################

my %Global_Aliases=();
my $ws="kbase";
my $bioObj = $FBAImpl->_get_msobject("Biochemistry",$ws,"default");
my $cpd_aliases = $bioObj->compoundsByAlias();

foreach my $aliasSet (sort keys %$cpd_aliases){
    foreach my $alias (sort keys %{$cpd_aliases->{$aliasSet}}){
	foreach my $id (sort keys %{$cpd_aliases->{$aliasSet}{$alias}}){
	    $Global_Aliases{$aliasSet}{$alias}{default}{$id}=1;
	}
    }
}

my $rxn_aliases = $bioObj->reactionsByAlias();

foreach my $aliasSet (sort keys %$rxn_aliases){
    foreach my $alias (sort keys %{$rxn_aliases->{$aliasSet}}){
	foreach my $id (sort keys %{$rxn_aliases->{$aliasSet}{$alias}}){
	    $Global_Aliases{$aliasSet}{$alias}{default}{$id}=1;
	}
    }
}

$bioObj = $FBAImpl->_get_msobject("Biochemistry",$ws,"plantdefault");
$cpd_aliases = $bioObj->compoundsByAlias();

#special case where previously merged compounds need to have their aliases restored
my %Merged_MS=();
foreach my $alias (keys %{$cpd_aliases->{ModelSEED}}){
    foreach my $id (keys %{$cpd_aliases->{ModelSEED}{$alias}}){
	if($alias ne $id){
	    $Merged_MS{$id}{$alias}=1;
	}
    }
}

foreach my $aliasSet (sort keys %$cpd_aliases){
    foreach my $alias (sort keys %{$cpd_aliases->{$aliasSet}}){
	foreach my $id (sort keys %{$cpd_aliases->{$aliasSet}{$alias}}){
	    $Global_Aliases{$aliasSet}{$alias}{plantdefault}{$id}=1;

	    if(exists($Merged_MS{$id})){
		foreach my $merged_id (keys %{$Merged_MS{$id}}){
		    $Global_Aliases{$aliasSet}{$alias}{plantdefault}{$merged_id}=1;
		}
	    }
	}
    }
}

$rxn_aliases = $bioObj->reactionsByAlias();

foreach my $aliasSet (sort keys %$rxn_aliases){
    foreach my $alias (sort keys %{$rxn_aliases->{$aliasSet}}){
	foreach my $id (sort keys %{$rxn_aliases->{$aliasSet}{$alias}}){
	    $Global_Aliases{$aliasSet}{$alias}{default}{$id}=1;
	}
    }
}

foreach my $aliasSet (sort keys %Global_Aliases){
    my $file = $aliasSet;
    $file = join("_",split(/\s/,$aliasSet)) if $aliasSet eq "Enzyme Class";

    open(OUT, "> ".$file.".aliases");
    print OUT $aliasSet."\tdefault\tplantdefault\n";
    foreach my $alias (sort keys %{$Global_Aliases{$aliasSet}}){
	print OUT $alias."\t";
	print OUT exists($Global_Aliases{$aliasSet}{$alias}{default}) ? join("|",sort keys %{$Global_Aliases{$aliasSet}{$alias}{default}}) : "";
	print OUT "\t";
	print OUT exists($Global_Aliases{$aliasSet}{$alias}{plantdefault}) ? join("|",sort keys %{$Global_Aliases{$aliasSet}{$alias}{plantdefault}}) : ""; 
	print OUT "\n";
    }
}
