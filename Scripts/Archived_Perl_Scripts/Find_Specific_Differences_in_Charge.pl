#!/usr/bin/env perl
use warnings;
use strict;
use Formulas;
my @temp=();

my $File = $ARGV[0];
die("No file or $File isn't found") if !$File || !-f $File;

#recover modifications that have already been applied to avoid duplicating records
#recover formerly modified formulas to automatically assign charges
open(FH, "< ../Biochemistry/compounds.master.mods");
my %Touched=();
my %Formulaed=();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $Touched{$temp[0]}=1 if $temp[2] eq "charge" ;
    $Formulaed{$temp[0]}=$temp[1];
}
close(FH);

open(FH, "< $File");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next if exists($Touched{$temp[0]});

    #remove single decimal place if zero
    $temp[1] =~ s/\.0$//;
    $temp[2] =~ s/\.0$//;

    #1 is default; 2 is plantdefault
    my $charge1 = $temp[1];
    my $charge2 = $temp[2];

    #If previously using formula, then adopt modifications
    #unless charge = 10000000
    if(exists($Formulaed{$temp[0]})){
	if($Formulaed{$temp[0]} eq "default" && $charge1 != 10000000){
	    print $temp[0],"\tdefault\tcharge\t$charge1\n";
	}elsif($Formulaed{$temp[0]} eq "plantdefault" && $charge2 != 10000000){
	    print $temp[0],"\tplantdefault\tcharge\t$charge2\n";
	}
	next;
    }

    #If missing charge
    if($charge1 == 10000000 || $charge2 == 10000000){
	if($charge1 == 10000000 && $charge2 != 10000000){
	    print $temp[0],"\tplantdefault\tcharge\t$charge2\n";
	}elsif($charge1 != 10000000 && $charge2 == 10000000){
	    print $temp[0],"\tdefault\tcharge\t$charge1\n";
	}
	next;
    }

#    print $temp[0],"\t",$charge1,"\t",$charge2,"\n";

}
