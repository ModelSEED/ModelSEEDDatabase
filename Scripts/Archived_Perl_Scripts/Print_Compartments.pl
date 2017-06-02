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
my @compartments = sort { $a->{id} cmp $b->{id} } @{$bioObj->compartments()};
open(OUT, "> ../Biochemistry/compartments.default.tsv");
print OUT "id\tname\thierarchy\n";
print OUT join("\n", map { $_->id()."\t".$_->name()."\t".$_->hierarchy() } @compartments),"\n";
close OUT;

$bioObj = $FBAImpl->_get_msobject("Biochemistry",$ws,"plantdefault");
@compartments = sort { $a->{id} cmp $b->{id} } @{$bioObj->compartments()};
open(OUT, "> ../Biochemistry/compartments.plantdefault.tsv");
print OUT "id\tname\thierarchy\n";
print OUT join("\n", map { $_->id()."\t".$_->name()."\t".$_->hierarchy() } @compartments),"\n";
close OUT;

$bioObj = $FBAImpl->_get_msobject("Biochemistry",$ws,"plantdefault_obs");
@compartments = sort { $a->{id} cmp $b->{id} } @{$bioObj->compartments()};
open(OUT, "> ../Biochemistry/compartments.plantdefault_obs.tsv");
print OUT "id\tname\thierarchy\n";
print OUT join("\n", map { $_->id()."\t".$_->name()."\t".$_->hierarchy() } @compartments),"\n";
close OUT;
