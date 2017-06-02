#!/bin/env perl
use strict;
use warnings;
my @temp=();

open(FH, "< Status_Diffs.txt");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);

    next if $temp[0] ne 1;

    my %Old_Codes = map { $_ => 1 } grep { length($_)>1 } map { (split(/:/,$_))[0] } split(/\|/,$temp[2]);
    my %New_Codes = map { $_ => 1 } grep { length($_)>1 } map { (split(/:/,$_))[0] } split(/\|/,$temp[3]);

#    print $_,"\n" if scalar(keys %Old_Codes)==0;
    my $Old = join("|", sort keys %Old_Codes);
    my $New = join("|", sort keys %New_Codes);

    print $temp[1],"\t",$Old,"\t",$New,"\n" if $Old ne $New;
}
