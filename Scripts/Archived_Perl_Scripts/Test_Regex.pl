#!/usr/bin/env perl
use warnings;
use strict;

my $String = '-1:cpd00001:c:0:"H2O";-1:cpd00102:c:0:"Glyceraldehyde3-phosphate";-2:cpd11621:c:0:"Oxidizedferredoxin";3:cpd00067:c:0:"H+";1:cpd00169:c:0:"3-Phosphoglycerate";2:cpd11620:c:0:"Reducedferredoxin"';

my $Cpd = 'cpd11620';
my $Stoich = '1';

$String =~ s/\d+:${Cpd}/${Stoich}:${Cpd}/;

print $String,"\n";
