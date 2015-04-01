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

my $bioObj = $FBAImpl->_get_msobject("Biochemistry","kbase","default");

my %Rxns_Codes=();
my %Codes_Rxns=();
my %Rxns=();
foreach my $rxn (@{$bioObj->reactions()}){
    next if $rxn->id() eq "rxn00000";

    my $code = $rxn->genEquationCode();
    if(!exists($Codes_Rxns{$code})){
	$code = $rxn->revGenEquationCode();
    }

    $Codes_Rxns{$code}{$rxn->id()}=1;
    $Rxns_Codes{$rxn->id()}{$code}=1;
    $Rxns{$rxn->id()}=$rxn;
}

#print join("\n", grep { scalar(keys %{$Rxns_Codes{$_}})>1 } sort keys %Rxns_Codes),"\n";

#foreach my $rxn (grep { scalar(keys %{$Rxns_Codes{$_}})>1 } sort keys %Rxns_Codes){
#    print $rxn,"\t",$Rxns{$rxn}->genEquation(),"\n";
#    foreach my $code (sort keys %{$Rxns_Codes{$rxn}}){
#	print "\t",$code,"\n";
#	foreach my $paired_rxn ( grep { $_ ne $rxn } sort keys %{$Codes_Rxns{$code}} ){
#	    print "\t\t",$paired_rxn,"\t",$Rxns{$paired_rxn}->genEquation(),"\n";
#	}
#    }
#    print "=================================================\n\n";
#}

#print join("\n",sort keys %Rxns_Codes),"\n";

my %Touched_Rxns=();
my %Touched_Codes=();
open(OUT, "> ../Biochemistry/reactions.default.cf.tsv");
print OUT "id\tprimary_name\tequation\tdefinition\tstatus\tpaired reactions\n";
foreach my $rxn (sort keys %Rxns_Codes){
    next if exists($Touched_Rxns{$rxn});

    my %Paired_Rxns=();
    foreach my $code (sort keys %{$Rxns_Codes{$rxn}}){
	$Touched_Rxns{$rxn}{$code}=1;
	$Touched_Codes{$code}{$rxn}=1;

	foreach my $paired_rxn (grep { $_ ne $rxn } sort keys %{$Codes_Rxns{$code}}){
	    $Paired_Rxns{$paired_rxn}=1;
	    $Touched_Rxns{$paired_rxn}{$code}=1;
	    $Touched_Codes{$code}{$paired_rxn}=1;
	}
    }

    print OUT $rxn."\t".$Rxns{$rxn}->name()."\t".$Rxns{$rxn}->genEquation()."\t".$Rxns{$rxn}->definition()."\t".$Rxns{$rxn}->status()."\t".join("|",sort keys %Paired_Rxns)."\n";
}
close(OUT);
