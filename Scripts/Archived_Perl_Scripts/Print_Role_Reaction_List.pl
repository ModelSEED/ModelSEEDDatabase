#!/usr/bin/env perl
use warnings;
use strict;
my $header=1;
my @temp=();

#KEGG Aliases
my %Aliases=();
open(FH, "< ../Aliases/KEGG.aliases");
$header=1;
while(<FH>){
    if($header){$header--;next}
    next unless $_ =~ /rxn/;
    chomp;
    @temp=split(/\t/,$_,-1);
    $Aliases{$temp[0]}= [ split(/\|/,$temp[1]) ];
}
close(FH);

#Pathways
my %Pathways=();
open(FH, "< ../Pathways/KEGG.pathways");
$header=1;
while(<FH>){
    if($header){$header--;next}
    chomp;
    @temp=split(/\t/,$_,-1);
    foreach my $kg_rxn (split(/\|/,$temp[4])){
	foreach my $ms_rxn (@{$Aliases{$kg_rxn}}){
	    $Pathways{$ms_rxn}{$temp[2]}=1;
	}
    }
}
close(FH);

#Complexes
my %Complexes=();
open(FH, "< ../Templates/Complexes.tsv");
$header=1;
while(<FH>){
    if($header){$header--;next}
    chomp;
    @temp=split(/\t/,$_,-1);
    foreach my $fr (split(/\|/,$temp[5])){
	$fr = ( split(/\;/,$fr) )[0];
	next if $fr eq "null";
	$Complexes{$temp[0]}{$fr}=1;
    }
}
close(FH);

#Roles
my %Roles=();
open(FH, "< ../Templates/Roles.tsv");
$header=1;
while(<FH>){
    if($header){$header--;next}
    chomp;
    @temp=split(/\t/,$_,-1);
    $Roles{$temp[0]}={name=>$temp[1]};
}
close(FH);

#Reactions
my %Reactions=();
open(FH, "< ../Templates/GramNegative/Reactions.tsv");
$header=1;
while(<FH>){
    if($header){$header--;next}
    chomp;
    @temp=split(/\t/,$_,-1);
    foreach my $cpx (split(/\|/,$temp[8])){
	foreach my $role (keys %{$Complexes{$cpx}}){
	    $Roles{$role}{reactions}{$temp[0]}=1;
	}
    }    
}
close(FH);

foreach my $role ( 
    grep { scalar(keys %{$Roles{$_}{reactions}})>0 }
    sort keys %Roles){
    foreach my $rxn (sort keys %{$Roles{$role}{reactions}}){
	print $Roles{$role}{name},"\t",$rxn,"\t",join("|",sort keys %{$Pathways{$rxn}}),"\n";
    }
}
