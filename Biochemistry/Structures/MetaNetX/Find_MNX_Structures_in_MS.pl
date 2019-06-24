#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();

open(FH, "< ../All_ModelSEED_Structures.txt");
my %MS_Structures=();
while(<FH>){
    chomp;
    @temp=split(/\t/);
    next unless $temp[1] eq "InChIKey";

    $temp[7] =~ s/-\w$//;
    $MS_Structures{$temp[7]}=1;
}
close(FH);

#Downloaded from:
#https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv
#06/21/2019
open(FH, "< chem_prop.tsv");
open(OUT, "> Structures_in_ModelSEED.txt");
while(<FH>){
    chomp;
    next if $_ =~ /^#/;
    @temp=split(/\t/,$_);
    next if !$temp[8];

    my $struct = $temp[8];
    $temp[8] =~ s/-\w$//;
    next if !exists($MS_Structures{$temp[8]});

    print(OUT $temp[0],"\t",$struct,"\n");
}
close(FH);

