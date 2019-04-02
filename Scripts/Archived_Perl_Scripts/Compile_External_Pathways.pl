#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my @temp=();

my $AliasRoot = "../../Biochemistry/Aliases/";
my $PwyRoot = $AliasRoot."Provenance/Primary_Databases/";

open(FH, "< ../../Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt");
my %Aliases_Reactions=();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $Aliases_Reactions{$temp[2]}{$temp[1]}{$temp[0]}=1;

    if($temp[2] eq "KEGG"){
	if($temp[1] =~ /_/){
	    $temp[1] = (split(/_/,$temp[1]))[0];
	    $Aliases_Reactions{$temp[2]}{$temp[1]}{$temp[0]}=1;
	}
    }

    if($temp[2] eq "MetaCyc"){
	if($temp[1] =~ /exp/){
	    $temp[1]=(split(/\.[a-z]+(\.metaexp)?/,$temp[1]))[0];
	    $Aliases_Reactions{$temp[2]}{$temp[1]}{$temp[0]}=1;
	}elsif($temp[1] =~ /\.[a-z]$/){
	    $temp[1] =~ s/\.[a-z]$//;
	    $Aliases_Reactions{$temp[2]}{$temp[1]}{$temp[0]}=1;
	}
    }
}
close(FH);

my %Lines=();
foreach my $biochem ("KEGG","MetaCyc"){
    open(FH, "< ".$PwyRoot."/".$biochem."_Pathways.tbl");
    while(<FH>){
	chomp;
	@temp=split(/\t/,$_,-1);
	if(exists($Aliases_Reactions{$biochem}) && exists($Aliases_Reactions{$biochem}{$temp[0]})){
	    foreach my $ms_rxn (keys %{$Aliases_Reactions{$biochem}{$temp[0]}}){
		my $line = $ms_rxn."\t".$temp[1]." (".$temp[2].")\t".$biochem;
		$Lines{$ms_rxn}{$line}=1;
	    }
	}
    }
}
open(OUT, "> ".$AliasRoot."Unique_ModelSEED_Reaction_Pathways.txt");
print OUT "ModelSEED ID\tExternal ID\tSource\n";
foreach my $rxn (sort keys %Lines){
    foreach my $line (sort keys %{$Lines{$rxn}}){
	print OUT $line."\n";
    }
}
