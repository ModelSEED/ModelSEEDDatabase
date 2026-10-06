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
close(OUT);

#Downloaded from:
#https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv
#06/21/2019
open(FH, "< chem_prop.tsv");
open(OUT, "> Structures_in_ModelSEED.txt");
my %MS_MNX_Structures=();
while(<FH>){
    chomp;
    next if $_ =~ /^#/;
    @temp=split(/\t/,$_);
    next if !$temp[8];

    my $struct = $temp[8];
    $temp[8] =~ s/-\w$//;
    next if !exists($MS_Structures{$temp[8]});

    $MS_MNX_Structures{$temp[0]}=$temp[8];
    print(OUT $temp[0],"\t",$struct,"\n");
}
close(FH);
close(OUT);

#Contrived from equilibrator's cache, see README.md
#06/21/2019
open(FH, "< eq_cpds.dump");
open(OUT, "> Structures_in_ModelSEED_and_eQuilibrator.txt");
while(<FH>){
    chomp;
    my ($mnx_id,$struct)=split(/\t/,$_);
    next if !$struct;

    my $temp=$struct;
    $temp =~ s/-\w$//;
    next if !exists($MS_MNX_Structures{$mnx_id}) || $MS_MNX_Structures{$mnx_id} ne $temp;

    print(OUT $mnx_id,"\t",$struct,"\n");
}
close(FH);
close(OUT);
