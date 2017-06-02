#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header = 1;

use Bio::KBase::ObjectAPI::KBaseOntology::Mapping;

open(DATA, "< ../Mappings/Mapping_Data.txt");
my @Data_Headers = split(/\t/,<DATA>,-1);
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

    if(-f $Map_Dir."Mapping_Subsystems.txt"){
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
    }

    $Mapping = Bio::KBase::ObjectAPI::KBaseOntology::Mapping->new($Mapping);
    open(OUT, "> ".$Map_Dir."Mapping.JSON");
    print OUT $Mapping->toJSON({pp=>1}),"\n";
    close(OUT);
}
close(DATA);
