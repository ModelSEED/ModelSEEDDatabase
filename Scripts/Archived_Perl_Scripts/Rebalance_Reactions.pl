#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my @temp=();
my $header=1;

use lib '../lib';
use Formulas;
use Reactions;

#Read status and stoichiometry
open(FH, "< ../Biochemistry/reactions.master.tsv");
my %Rxn_Attrs=();
$header=1;
while(<FH>){
    chomp;
    if($header){$header--;next}
    @temp=split(/\t/,$_,-1);
    $Rxn_Attrs{$temp[0]}{"id"}=$temp[0];
    $Rxn_Attrs{$temp[0]}{"status"}=$temp[17];
    $Rxn_Attrs{$temp[0]}{"stoichiometry"}=$temp[4];
    $Rxn_Attrs{$temp[0]}{"direction"}=$temp[8];
    $Rxn_Attrs{$temp[0]}{"direction"}=$temp[9] if $temp[8] eq "?";
}
close(FH);

#Read name,charge, formula
open(FH, "< ../Biochemistry/compounds.master.tsv");
my %Cpd_Attrs=();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $Cpd_Attrs{$temp[0]}{"id"}=$temp[0];
    $Cpd_Attrs{$temp[0]}{"name"}=$temp[2];
    $Cpd_Attrs{$temp[0]}{"formula"}=$temp[3];
    $Cpd_Attrs{$temp[0]}{"charge"}=$temp[7];
}
close(FH);

my %Rxns_Cpds=();
my %Cpds_Rxns=();
foreach my $rxn (keys %Rxn_Attrs){
    foreach my $rgt (split(/;/,$Rxn_Attrs{$rxn}{'stoichiometry'})){
	my ($coeff,$rgt_id,$cpt_id,$index,$name) = split(/:/,$rgt);
	$Rxns_Cpds{$rxn}{$rgt_id}={"coefficient"=>$coeff,"compartment"=>$cpt_id};
	$Cpds_Rxns{$rgt_id}{$rxn}={"coefficient"=>$coeff,"compartment"=>$cpt_id};
    }
}

foreach my $rxn (sort keys %Rxn_Attrs){
    next unless $Rxn_Attrs{$rxn}->{status} eq "TOBAL";

    #Build Cpd Hash
    my $Contains_Null=0;
    my $Cpd_Hash=();
    foreach my $cpd (sort keys %{$Rxns_Cpds{$rxn}}){
	if($Cpd_Attrs{$cpd}{'formula'} ne "null"){
	    print $cpd,"\t",$Cpd_Attrs{$cpd}{'formula'},"\t",$Rxns_Cpds{$rxn}{$cpd}{'coefficient'},"\t",$Cpd_Attrs{$cpd}{'charge'},"\n";
	    $Cpd_Hash->{$cpd}={'formula'=> $Cpd_Attrs{$cpd}{'formula'},
			       'charge' => $Cpd_Attrs{$cpd}{'charge'},
			       'coefficient' => $Rxns_Cpds{$rxn}{$cpd}{'coefficient'},
			       'compartment' => $Rxns_Cpds{$rxn}{$cpd}{'compartment'}};
	}else{
	    $Contains_Null=1;
	}
    }

    my $Status = Reactions::balance_reaction($Cpd_Hash);
    print $Status,"\n";
}
