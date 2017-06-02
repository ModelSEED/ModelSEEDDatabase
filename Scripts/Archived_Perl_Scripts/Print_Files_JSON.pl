#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header = 1;

open(FH, "< ../Biochemistry/reactions.master.tsv");
my @headers = split(/\t/,<FH>);
chomp($headers[$#headers]);
my %Rxns=();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    my %Rxn_Hash = map { $headers[$_] => $temp[$_] } (0..$#temp);
    $Rxns{$temp[0]}=\%Rxn_Hash;
}
close(FH);

use JSON;
open(OUT, "> ../Biochemistry/reactions.master.json");
print OUT to_json(\%Rxns,{pretty=>1});
close(OUT);
