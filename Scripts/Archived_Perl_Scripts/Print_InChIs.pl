#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my @temp=();
my $BPath = "../../Biochemistry/";
my %Structures=();

open(FH, "< ".$BPath."Structures/KEGG_Original_InChI.txt");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_);
    $Structures{$temp[0]}=$temp[1];
}
close(FH);

open(FH, "< ".$BPath."Structures/MetaCyc_Original_InChI.txt");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_);
    $Structures{$temp[0]}=$temp[1];
}
close(FH);

my %Global_Structures=();
open(FH, "< ".$BPath."Aliases/Compounds_Aliases.tsv");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next unless $temp[3] =~ /KEGG|MetaCyc/;
    next unless exists($Structures{$temp[2]});

    foreach my $cpd (split(/\|/,$temp[0])){
	$Global_Structures{$temp[3]}{$temp[2]}{$Structures{$temp[2]}}{$cpd}=1;
    }
    if($temp[1]){
	$Global_Structures{$temp[3]}{$temp[2]}{$Structures{$temp[2]}}{$temp[1]}=1;
    }
}

my $file = "Structures.tsv";
open(MAIN, "> ".$file);
print MAIN "MS ID\tExternal ID\tSource\tStructure\n";

foreach my $aliasSet (sort keys %Global_Structures){
    foreach my $alias (sort keys %{$Global_Structures{$aliasSet}}){
	foreach my $structure (sort keys %{$Global_Structures{$aliasSet}{$alias}}){
	    print MAIN join("|",sort keys %{$Global_Structures{$aliasSet}{$alias}{$structure}}),"\t";
	    print MAIN $alias,"\t",$aliasSet,"\t",$structure,"\n";
	}
    }
}
close(MAIN);
