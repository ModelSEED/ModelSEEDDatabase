#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my @temp=();
my $BPath = "../../Biochemistry/";
my %Structures=();

open(FH, "< ".$BPath."Structures/KEGG/InChI_ChargedStrings.txt");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_);
    $Structures{$temp[0]}{'InChI'}=$temp[1];
}
close(FH);

open(FH, "< ".$BPath."Structures/MetaCyc/InChI_ChargedStrings.txt");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_);
    $Structures{$temp[0]}{'InChI'}=$temp[1];
}
close(FH);

open(FH, "< ".$BPath."Structures/KEGG/SMILE_ChargedStrings.txt");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_);
    $Structures{$temp[0]}{'SMILE'}=$temp[1];
}
close(FH);

open(FH, "< ".$BPath."Structures/MetaCyc/SMILE_ChargedStrings.txt");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_);
    $Structures{$temp[0]}{'SMILE'}=$temp[1];
}
close(FH);

my %Global_Structures=();
open(FH, "< ".$BPath."Aliases/Compounds_Aliases.tsv");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next unless $temp[3] =~ /KEGG|MetaCyc/;
    next unless exists($Structures{$temp[2]}) && exists($Structures{$temp[2]}{'InChI'});

    foreach my $cpd (split(/\|/,$temp[0])){
	$Global_Structures{$temp[3]}{$temp[2]}{$Structures{$temp[2]}{'InChI'}}{$cpd}=1;
    }
    if($temp[1]){
	$Global_Structures{$temp[3]}{$temp[2]}{$Structures{$temp[2]}{'InChI'}}{$temp[1]}=1;
    }
}

my $file = $BPath."/Structures/Structures.tsv";
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
