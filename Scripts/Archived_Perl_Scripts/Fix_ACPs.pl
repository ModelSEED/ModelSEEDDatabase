#!/usr/bin/env perl
use warnings;
use strict;
use Formulas;

my @temp=();
my @Lines=();
my $PP_Formula = {C=>11,H=>21,N=>2,O=>7,P=>1}; #Taken from Charged InChI string for C01134
my $PP_Charge = -2;

open(FH, "< Prioritized_ACPs_3.txt");
open(OUT, "> Prioritized_ACPs_3.fixes");
#open(FH, "< Prioritized_acyl-carriers.txt");
#open(OUT, "> Prioritized_acyl-carriers.fixes");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);

    my $formula = Formulas::parse($temp[2]);

    #All ACPs updated to charge of PP group
    print OUT $temp[0],"\tcurated\tcharge\t-2\n";

    next if exists($formula->{'P'});
    
    foreach my $atom (keys %$PP_Formula){
	$formula->{$atom} += $PP_Formula->{$atom};
    } 

    print $temp[0],"\t",$temp[2],"\t",Formulas::hill_sort_formula($formula),"\n";

    print OUT $temp[0],"\tcurated\tformula\t".Formulas::hill_sort_formula($formula)."\n";

}
