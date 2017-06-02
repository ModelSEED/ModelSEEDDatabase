#!/usr/bin/env perl
use warnings;
use strict;

#######################################################
#Initialization
#######################################################

my $directory = $ARGV[0];

use Bio::KBase::fbaModelServices::ScriptHelpers qw( getToken );
my $AToken = getToken();

use Bio::KBase::fbaModelServices::Impl;
my $fba = Bio::KBase::fbaModelServices::Impl->new({'fbajobcache' => "/homes/seaver/Projects/KBase_Scripts/FBA_Scripts/JobsCache",
						       'jobserver-url' => "http://kbase.us/services/workspace",
						       'fbajobdir' => "/tmp/fbajobs",
						       'mfatoolkitbin' => "/vol/model-prod/kbase/MFAToolkit/bin/mfatoolkit",
						       'probanno-url' => "http://kbase.us/services/probabilistic_annotation/",
						       'mssserver-url' => "http://bio-data-1.mcs.anl.gov/services/ms_fba",
						       'accounttype' => "kbase",
						       'workspace-url' => "http://kbase.us/services/ws",
						       'defaultJobState' => "queued",
						       'gaserver-url' => "http://kbase.us/services/genome_annotation",
						       'idserver-url' => "http://kbase.us/services/idserver"});
$fba->_setContext(undef,{auth=>$AToken});

# Get the default Mapping object.
my $map = $fba->_get_msobject("Mapping","kbase","default-mapping");

#Printing complex roles
open(my $fh, ">", $directory."/ComplexRoles.default.tsv");
my $columns = [qw(
	complex_id
	complex_name
	complex_source
	complex_type
	role_id
	role_name
	role_type
	role_source
	role_aliases
	role_exemplar
	type
	triggering
	optional
)];
print $fh join("\t",@{$columns})."\n";
my $cpxs = $map->complexes();
my $idhash;
for (my $i=0; $i < @{$cpxs}; $i++) {
	my $cpx = $cpxs->[$i];
	$idhash->{$cpx->id()} = "cpx.".$i;
	my $cpxroles = $cpx->complexroles();
	for (my $j=0; $j < @{$cpxroles}; $j++) {
		my $cpxrole = $cpxroles->[$j];
		my $roleid = $cpxrole->role()->id();
		$roleid =~ s/ms//;
		my $data = [
	    	"cpx.".$i,
	    	"cpx.".$i,
	    	"ModelSEED",
	    	"SEED_role_complex",
	    	$roleid,
	    	$cpxrole->role()->name(),
	    	"SEED_role",
	    	"SEED",
	    	"searchname:".$cpxrole->role()->searchname(),
	    	"null",
	    	"role_mapping",
	    	$cpxrole->triggering(),
	    	$cpxrole->optionalRole(),
	    ];
	    print $fh join("\t",@{$data})."\n";
	}
}
close($fh);

