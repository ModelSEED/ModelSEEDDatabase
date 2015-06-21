#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header = 1;

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

open(DATA, "< ../ModelTemplates/ModelTemplate_Data.txt");
my @Data_Headers = split(/\t/,<DATA>,-1);
chomp($Data_Headers[$#Data_Headers]);
while(<DATA>){
    chomp;
    @temp=split(/\t/,$_,-1);
    my $MT_Dir = "../ModelTemplates/".$temp[0]."/";
    my $ModelTemplate = {};

    for(my $i=1;$i<scalar(@Data_Headers);$i++){
	$ModelTemplate->{$Data_Headers[$i]}=$temp[$i];
    }

    open(TMPLRXN, "< ".$MT_Dir."ModelTemplate_Reactions.txt");
    my @TmplRxn_Headers = split(/\t/,<TMPLRXN>,-1);
    chomp($TmplRxn_Headers[$#TmplRxn_Headers]);
    while(<TMPLRXN>){
	chomp;
	@temp = split(/\t/,$_,-1);
	my $TmplRxn = {};
	for(my $i=0;$i<scalar(@TmplRxn_Headers);$i++){
	    if($TmplRxn_Headers[$i] eq "complex_refs"){
		$TmplRxn->{$TmplRxn_Headers[$i]}=[split(/\|/,$temp[$i])];
	    }else{
		$TmplRxn->{$TmplRxn_Headers[$i]}=$temp[$i];
	    }
	}
	push(@{$ModelTemplate->{templateReactions}},$TmplRxn);
    }
    close(TMPLRXN);

    open(TMPLBIO, "< ".$MT_Dir."ModelTemplate_Biomasses.txt");
    my @TmplBio_Headers = split(/\t/,<TMPLBIO>,-1);
    chomp($TmplBio_Headers[$#TmplBio_Headers]);
    while(<TMPLBIO>){
	chomp;
	@temp = split(/\t/,$_,-1);
	my $TmplBio = {};
	for(my $i=0;$i<scalar(@TmplBio_Headers);$i++){
	    $TmplBio->{$TmplBio_Headers[$i]}=$temp[$i];
	}
	open(TMPLBIOCMP, "< ".$MT_Dir."ModelTemplate_Biomass_Components.txt");
	my @TmplBioCmp_Headers = split(/\t/,<TMPLBIOCMP>,-1);
	chomp($TmplBioCmp_Headers[$#TmplBioCmp_Headers]);
	while(<TMPLBIOCMP>){
	    chomp;
	    @temp = split(/\t/,$_,-1);

	    #Only biomass components whose id matches
	    my $TmplBio_ID = $TmplBio->{id};
	    next unless $temp[0] =~ /^${TmplBio_ID}/;

	    my $TmplBioCmp = {};
	    for(my $i=0;$i<scalar(@TmplBioCmp_Headers);$i++){
		if($TmplBioCmp_Headers[$i] eq "linked_compound_refs" || $TmplBioCmp_Headers[$i] eq "link_coefficients"){
		    $TmplBioCmp->{$TmplBioCmp_Headers[$i]}=[split(/\|/,$temp[$i])];
		}else{
		    $TmplBioCmp->{$TmplBioCmp_Headers[$i]}=$temp[$i];
		}
	    }
	    push(@{$TmplBio->{templateBiomassComponents}},$TmplBioCmp);
	}
	push(@{$ModelTemplate->{templateBiomasses}},$TmplBio);
    }
    close(TMPLBIO);

    $ModelTemplate = Bio::KBase::ObjectAPI::KBaseFBA::ModelTemplate->new($ModelTemplate);
    open(OUT, "> ".$MT_Dir."ModelTemplate.JSON");
    print OUT $ModelTemplate->toJSON({pp=>1}),"\n";
    close(OUT);
}
close(DATA);

open(DATA, "< ../Mappings/Mapping_Data.txt");
@Data_Headers = split(/\t/,<DATA>,-1);
chomp($Data_Headers[$#Data_Headers]);
while(<DATA>){
    chomp;
    @temp=split(/\t/,$_,-1);
    my $Map_Dir = "../Mappings/".$temp[0]."/";
    my $Mapping = {};

    print $Map_Dir,"\n";
    for(my $i=1;$i<scalar(@Data_Headers);$i++){
	if($Data_Headers[$i] =~ /aliases$/){
	    next if !$temp[$i];

	    my $hash_ref = {};
	    foreach my $key (split(/\|/,$temp[$i])){
		my ($id,$alias_sets) = split(/:/,$key,2);
		foreach my $alias_set (split(/\//,$alias_sets)){
		    my ($set_id,$aliases) = split(/:/,$alias_set);
		    $hash_ref->{$id}{$set_id}=[split(/;/,$aliases)];
		}
	    }
	    $Mapping->{$Data_Headers[$i]}=$hash_ref;
	}else{
	    $Mapping->{$Data_Headers[$i]}=$temp[$i];
	}
    }

    open(ROLES, "< ".$Map_Dir."Mapping_Roles.txt");
    my @Role_Headers = split(/\t/,<ROLES>,-1);
    chomp($Role_Headers[$#Role_Headers]);
    while(<ROLES>){
	chomp;
	@temp=split(/\t/,$_,-1);

	my $Role = {};
	for(my $i=1;$i<scalar(@Role_Headers);$i++){
	    $Role->{$Role_Headers[$i]}=$temp[$i];
	}

	push(@{$Mapping->{roles}},$Role);
    }
    close(ROLES);

    open(SUBSYSTEMS, "< ".$Map_Dir."Mapping_Subsystems.txt");
    my @Subsystem_Headers = split(/\t/,<SUBSYSTEMS>,-1);
    chomp($Subsystem_Headers[$#Subsystem_Headers]);
    while(<SUBSYSTEMS>){
	chomp;
	@temp=split(/\t/,$_,-1);

	my $Subsystem = {};
	for(my $i=1;$i<scalar(@Subsystem_Headers);$i++){
	    if($Subsystem_Headers[$i] eq "role_refs"){
		$Subsystem->{$Subsystem_Headers[$i]}=[split(/\|/,$temp[$i])];
	    }else{
		$Subsystem->{$Subsystem_Headers[$i]}=$temp[$i];
	    }
	}

	push(@{$Mapping->{subsystems}},$Subsystem);
    }
    close(SUBSYSTEMS);

    open(COMPLEXES, "< ".$Map_Dir."Mapping_Complexes.txt");
    my @Complex_Headers = split(/\t/,<COMPLEXES>,-1);
    chomp($Complex_Headers[$#Complex_Headers]);
    while(<COMPLEXES>){
	chomp;
	@temp=split(/\t/,$_,-1);

	my $Complex = {};
	for(my $i=1;$i<scalar(@Complex_Headers);$i++){
	    if($Complex_Headers[$i] eq "complexroles"){
		my $hash_ref = { map { my @Pair = split(/:/,$_); $Pair[0] => $Pair[1] } split(/;/,$temp[$i]) };
		push(@{$Complex->{$Complex_Headers[$i]}},$hash_ref);
	    }else{
		$Complex->{$Complex_Headers[$i]}=$temp[$i];
	    }
	}

	push(@{$Mapping->{complexes}},$Complex);
    }
    close(COMPLEXES);

    $Mapping = Bio::KBase::ObjectAPI::KBaseOntology::Mapping->new($Mapping);
    open(OUT, "> ".$Map_Dir."Mapping.JSON");
    print OUT $Mapping->toJSON({pp=>1}),"\n";
    close(OUT);
}
close(DATA);
