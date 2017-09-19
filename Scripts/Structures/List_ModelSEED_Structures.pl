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

my %Global_InChI_Structures=();
my %Global_SMILE_Structures=();
open(FH, "< ".$BPath."Aliases/Compounds_Aliases.tsv");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next unless $temp[3] =~ /KEGG|MetaCyc/;
    next unless exists($Structures{$temp[2]});

    foreach my $cpd (split(/\|/,$temp[0])){
	if(exists($Structures{$temp[2]}{'InChI'})){
	    $Global_InChI_Structures{$temp[3]}{$temp[2]}{$Structures{$temp[2]}{'InChI'}}{$cpd}=1;
	}
	if(exists($Structures{$temp[2]}{'SMILE'})){
	    $Global_SMILE_Structures{$temp[3]}{$temp[2]}{$Structures{$temp[2]}{'SMILE'}}{$cpd}=1;
	}
    }
    if($temp[1]){
	if(exists($Structures{$temp[2]}{'InChI'})){
	    $Global_InChI_Structures{$temp[3]}{$temp[2]}{$Structures{$temp[2]}{'InChI'}}{$temp[1]}=1;
	}
	if(exists($Structures{$temp[2]}{'SMILE'})){
	    $Global_SMILE_Structures{$temp[3]}{$temp[2]}{$Structures{$temp[2]}{'SMILE'}}{$temp[1]}=1;
	}
    }
}

my $file = $BPath."/Structures/Structures.tsv";
open(MAIN, "> ".$file);
print MAIN "MS ID\tExternal ID\tSource\tStructure\n";
foreach my $aliasSet (sort keys %Global_InChI_Structures){
    foreach my $alias (sort keys %{$Global_InChI_Structures{$aliasSet}}){
	foreach my $structure (sort keys %{$Global_InChI_Structures{$aliasSet}{$alias}}){
	    print MAIN join("|",sort keys %{$Global_InChI_Structures{$aliasSet}{$alias}{$structure}}),"\t";
	    print MAIN $alias,"\t",$aliasSet,"\t",$structure,"\n";
	}
    }
}
foreach my $aliasSet (sort keys %Global_SMILE_Structures){
    foreach my $alias (sort keys %{$Global_SMILE_Structures{$aliasSet}}){
	foreach my $structure (sort keys %{$Global_SMILE_Structures{$aliasSet}{$alias}}){
	    print MAIN join("|",sort keys %{$Global_SMILE_Structures{$aliasSet}{$alias}{$structure}}),"\t";
	    print MAIN $alias,"\t",$aliasSet,"\t",$structure,"\n";
	}
    }
}
close(MAIN);
