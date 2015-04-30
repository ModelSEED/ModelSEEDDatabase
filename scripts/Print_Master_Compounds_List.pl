#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header=1;

#Start with original biochemistry
open(FH, "< ../Biochemistry/compounds.default.tsv");
$header=1;
my %Cpds=();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $Cpds{$temp[0]}=$_;
}
close(FH);

#Add new biochemistry
open(FH, "< ../Biochemistry/compounds.plantdefault.tsv");
$header=1;
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $Cpds{$temp[0]}=$_ if !exists($Cpds{$temp[0]});
}
close(FH);

#Add aliases for KEGG and MetaCyc to match InChI structures
open(FH, "< ../Aliases/KEGG.aliases");
my %Aliases = ();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    foreach my $id (split(/\|/,$temp[2])){
	$Aliases{$id}{KEGG}=$temp[0];
    }
}
close(FH);

open(FH, "< ../Aliases/MetaCyc.aliases");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    foreach my $id (split(/\|/,$temp[2])){
	$Aliases{$id}{MetaCyc}=$temp[0];
    }
}
close(FH);

#Collect structures
open(FH, "< ../Structures/KEGG_Search_InChI.txt");
my %InChIs = ();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $InChIs{$temp[0]}=$temp[1];
}
close(FH);
open(FH, "< ../Structures/MetaCyc_Search_InChI.txt");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $InChIs{$temp[0]}=$temp[1];
}
close(FH);

#Print it all out
open(OUT, "> Master_Compound_List.tsv");
foreach my $cpd ( grep { $_ ne "cpd00000" } sort keys %Cpds){
    print OUT $Cpds{$cpd}."\t";

    #priortize InChIs from KEGG
    if(exists($Aliases{$cpd}{KEGG}) && exists($InChIs{$Aliases{$cpd}{KEGG}})){
	print OUT $InChIs{$Aliases{$cpd}{KEGG}}."\n";
    }elsif(exists($Aliases{$cpd}{MetaCyc}) && exists($InChIs{$Aliases{$cpd}{MetaCyc}})){
	print OUT $InChIs{$Aliases{$cpd}{MetaCyc}}."\n";
    }
}
