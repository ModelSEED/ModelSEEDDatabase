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

my $ws="kbase";
my $bioObj = $FBAImpl->_get_msobject("Biochemistry",$ws,"default");
open(OUT, "> ../Biochemistry/compounds.default.tsv");
print OUT "id\tname\tabbreviation\tformula\tcharge\tisCofactor\n";
print OUT join("\n", map { $_->id()."\t".$_->name()."\t".$_->abbreviation()."\t".$_->formula()."\t".$_->defaultCharge()."\t".$_->isCofactor() } @{$bioObj->compounds()}),"\n";
close OUT;

$bioObj = $FBAImpl->_get_msobject("Biochemistry",$ws,"plantdefault");
open(OUT, "> ../Biochemistry/compounds.plantdefault.tsv");
print OUT "id\tname\tabbreviation\tformula\tcharge\tisCofactor\n";
print OUT join("\n", map { $_->id()."\t".$_->name()."\t".$_->abbreviation()."\t".$_->formula()."\t".$_->defaultCharge()."\t".$_->isCofactor() } @{$bioObj->compounds()}),"\n";
close OUT;

$bioObj = $FBAImpl->_get_msobject("Biochemistry",$ws,"plantdefault_obs");
open(OUT, "> ../Biochemistry/compounds.plantdefault_obs.tsv");
print OUT "id\tname\tabbreviation\tformula\tcharge\tisCofactor\n";
print OUT join("\n", map { $_->id()."\t".$_->name()."\t".$_->abbreviation()."\t".$_->formula()."\t".$_->defaultCharge()."\t".$_->isCofactor() } @{$bioObj->compounds()}),"\n";
close OUT;
