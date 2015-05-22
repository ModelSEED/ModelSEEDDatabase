#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header=1;

#Load up the required modifications
open(FH, "< ../Biochemistry/reactions.master.mods");
my %Rxn_Mods=();
while(<FH>){
    chomp;
    @temp = split(/\t/,$_,-1);
    $Rxn_Mods{$temp[0]}{$temp[2]}=$temp[3];
}
close(FH);

#Start with original biochemistry
open(FH, "< ../Biochemistry/reactions.default.tsv");
my %Rxns=();
my @headers = split(/\t/,<FH>);
chomp $headers[$#headers];
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    for(my $i=1;$i<scalar(@temp);$i++){

	if(exists($Rxn_Mods{$temp[0]}) && exists($Rxn_Mods{$temp[0]}{$headers[$i]})){
	    $temp[$i] = $Rxn_Mods{$temp[0]}{$headers[$i]};
	}

	$Rxns{$temp[0]}{$headers[$i]}=$temp[$i];
    }
}
close(FH);

#Add new biochemistry
open(FH, "< ../Biochemistry/reactions.plantdefault.tsv");
$header=1;
while(<FH>){
    chomp;
    if($header){$header--;next;}
    @temp=split(/\t/,$_,-1);
    next if exists($Rxns{$temp[0]});

    for(my $i=1;$i<scalar(@temp);$i++){
	if(exists($Rxn_Mods{$temp[0]}) && exists($Rxn_Mods{$temp[0]}{$headers[$i]})){
	    $temp[$i] = $Rxn_Mods{$temp[0]}{$headers[$i]};
	}

	$Rxns{$temp[0]}{$headers[$i]}=$temp[$i];
    }
}
close(FH);

#Print it all out
open(OUT, "> Master_Reaction_List.tsv");
print OUT join("\t",@headers),"\n",;
foreach my $cpd ( grep { $_ ne "cpd00000" } sort keys %Rxns){
    print OUT $cpd."\t".join("\t", map { $Rxns{$cpd}{$_} } grep { $_ ne "id" } @headers),"\n";
}
