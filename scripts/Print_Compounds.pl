#!/usr/bin/env perl
use warnings;
use strict;

my $DB = $ARGV[0];
exit if !$DB;
exit if $DB ne "default" && $DB ne "plantdefault";

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
my $bioObj = $FBAImpl->_get_msobject("Biochemistry","kbase", $DB);

#List of headers used in SOLR Dump
my @headers = ("id","abbreviation","name","formula","mass","source","structure","charge","is_core","is_obsolete","linked_compound","is_cofactor","deltag","deltagerr","pka","pkb","abstract_compound","comprised_of","aliases");

#Translation of header to originating attribute in compound
my %Original_Header = ( charge=>"defaultCharge", deltag=>"deltaG", deltagerr => "deltaGErr", is_cofactor => "isCofactor" );

#Headers to skip as they are not true compound attributes in current biochemistry
my %Skip_Header = ( source => 1, structure => 1, is_core => 1, is_obsolete => 1, linked_compound => 1, comprised_of => 1, aliases => 1 );

#Default values to use
my %Default_Values = ( mass => "null", deltag => "null", deltagerr => "null", comprised_of => "null", source => "ModelSEED", is_core => 1, is_obsolete => 0,
		       linked_compound => "null", comprised_of => "null", structure => "null", aliases => "null" );

#Hash of subroutines that transform the values of some hierarchical attributes to a single string
#The subroutines here are copied from the DumpSOLRTables scripts
my %Transform_Header =  map { my $item = pop @$_; map { $_, $item } @$_ }
[qw(abstract_compound) => sub { my $cpd = shift; return defined($cpd->abstractCompound_ref()) ? $cpd->abstractCompound()->id() : "null"; }],
[qw(pka) => sub { my $cpd = shift; my $pka= join(";", map { my $atom = $_; map { $_.":".$atom } @{$cpd->pkas()->{$atom}} } keys %{$cpd->pkas()}); return $pka; }],
[qw(pkb) => sub { my $cpd = shift; my $pkb= join(";", map { my $atom = $_; map { $_.":".$atom } @{$cpd->pkbs()->{$atom}} } keys %{$cpd->pkbs()}); return $pkb; }];

my @compounds = sort { $a->{id} cmp $b->{id} } grep { $_->id() ne "cpd00000" } @{$bioObj->compounds()};

open(OUT, "> ../Biochemistry/compounds.".$DB.".tsv");
print OUT join("\t", @headers),"\n";
foreach my $cpd (@compounds){

#Here, we either transform the attributes value, use its original value (or skip it), and then use the default value if availeble
#Finally, the header itself if printed if it falls through to the end, so a simple test is to see if any column headers are repeated
#Lots of nulls appear at this stage, as they are integrated when forming the master document
print OUT join("\t", map { 
    if(exists($Transform_Header{$_})){ 
	$Transform_Header{$_}($cpd);
    }elsif( !exists($Skip_Header{$_}) && defined($cpd->$_()) ){
	$cpd->$_();
    }elsif( exists($Default_Values{$_}) ){ 
	$Default_Values{$_}; 
    }else{
	$_;
    }} map { exists($Original_Header{$_}) ? $Original_Header{$_} : $_ } @headers),"\n";
}
close(OUT);
