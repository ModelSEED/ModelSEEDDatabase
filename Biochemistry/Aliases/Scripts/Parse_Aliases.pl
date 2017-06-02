#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my @temp=();

opendir(my $dh, "../");
my @Files = grep { $_ =~ /\.aliases/ } readdir($dh);
closedir($dh);

my %Cpd_Mismatch=();
my %Cpd_Mismatch_Count=();

my %Rxn_Mismatch=();
my %Rxn_Mismatch_Count=();
foreach my $file (@Files){
    my $filestub = $file;
    $filestub =~ s/\.aliases//;

    open(FH, "< ../$file");
    while(<FH>){
	chomp;
	@temp=split(/\t/,$_,3);
	
	if($temp[1] =~ /cpd/ && $temp[2] =~ /cpd/ && $temp[1] ne $temp[2]){
	    foreach my $ms_cpd (split(/\|/,$temp[1])){
		foreach my $ps_cpd (split(/\|/,$temp[2])){
		    next if $ms_cpd eq $ps_cpd;
		    my $pair = $ms_cpd.":".$ps_cpd;
		    $Cpd_Mismatch{$pair}{$temp[0]}{$filestub}=1;
		}
	    }
	    $Cpd_Mismatch_Count{$filestub}++;
	}

	if($temp[1] =~ /rxn/ && $temp[2] =~ /rxn/ && $temp[1] ne $temp[2]){
	    foreach my $ms_rxn (split(/\|/,$temp[1])){
		foreach my $ps_rxn (split(/\|/,$temp[2])){
		    next if $ms_rxn eq $ps_rxn;
		    my $pair = $ms_rxn.":".$ps_rxn;
		    $Rxn_Mismatch{$pair}{$temp[0]}{$filestub}=1;
		}
	    }
	    $Rxn_Mismatch_Count{$filestub}++;
	}

    }
}

open(OUT, "> Cpd_Mismatch.json");
print OUT to_json(\%Cpd_Mismatch,{pretty=>1}),"\n";
close(OUT);

open(OUT, "> Cpd_Mismatch_Count.txt");
print OUT join("\n", map { $_."\t".$Cpd_Mismatch_Count{$_} } sort { $Cpd_Mismatch_Count{$b} <=> $Cpd_Mismatch_Count{$a} } keys %Cpd_Mismatch_Count),"\n";
close(OUT);

open(OUT, "> Rxn_Mismatch.json");
print OUT to_json(\%Rxn_Mismatch,{pretty=>1}),"\n";
close(OUT);

open(OUT, "> Rxn_Mismatch_Count.txt");
print OUT join("\n", map { $_."\t".$Rxn_Mismatch_Count{$_} } sort { $Rxn_Mismatch_Count{$b} <=> $Rxn_Mismatch_Count{$a} } keys %Rxn_Mismatch_Count),"\n";
close(OUT);
