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
my @headers = ("id","abbreviation","name","code","stoichiometry","is_transport","equation","definition","reversibility","direction","abstract_reaction","pathways","aliases","ec_numbers","deltag","deltagerr","compound_ids","status");

#Translation of header to originating attribute in reaction
my %Original_Header = ( deltag=>"deltaG", deltagerr => "deltaGErr", is_transport => "isTransport", reversibility => "thermoReversibility" );

#Headers to skip as they are not true compound attributes in current biochemistry
my %Skip_Header = ( abstract_reaction => 1, pathways => 1, aliases => 1, ec_numbers => 1, compound_ids => 1 );

#Default values to use
my %Default_Values = ( stoichiometry => "null", deltaG => "null", deltaGErr => "null", abstract_reaction => "null",
		       compound_ids => "null", aliases => "null", ec_numbers => "null", aliases => "null", pathways => "null" );

#Hash of subroutines that transform the values of some hierarchical attributes to a single string
#The subroutines here are copied from the DumpSOLRTables scripts
my %Transform_Header =  map { my $item = pop @$_; map { $_, $item } @$_ }
[qw(abstract_reaction) => sub { my $rxn = shift; return defined($rxn->abstractReaction_ref()) ? $rxn->abstractReaction()->id() : "null"; }],
[qw(compound_ids) => sub { my $rxn = shift; my $compounds = join(";", map { $_->compound()->id() } @{$rxn->reagents()}); return $compounds; }];

my @reactions = sort { $a->{id} cmp $b->{id} } grep { $_->id() ne "rxn00000" && $_->id() !~ /^CoA-2-methylpropanoylating/ } @{$bioObj->reactions()};

open(OUT, "> ../Biochemistry/reactions.".$DB.".tsv");
print OUT join("\t", @headers),"\n";
foreach my $rxn (@reactions){

#Here, we either transform the attributes value, use its original value (or skip it), and then use the default value if availeble
#Finally, the header itself if printed if it falls through to the end, so a simple test is to see if any column headers are repeated
#Lots of nulls appear at this stage, as they are integrated when forming the master document
print OUT join("\t", map { 
    if(exists($Transform_Header{$_})){ 
	$Transform_Header{$_}($rxn);
    }elsif( !exists($Skip_Header{$_}) && defined($rxn->$_()) ){
	$rxn->$_();
    }elsif( exists($Default_Values{$_}) ){ 
	$Default_Values{$_}; 
    }else{
	$_;
    }} map { exists($Original_Header{$_}) ? $Original_Header{$_} : $_ } @headers),"\n";
}
close(OUT);
