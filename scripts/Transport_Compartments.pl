#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header=1;

my %Main_Cmpts=();
open(FH, "< ../Biochemistry/reactions.plantdefault_obs.tsv");
$header=1;
while(<FH>){
    chomp;
    if($header){$header--;next}
    @temp=split(/\t/,$_,-1);

    my %Cmpts=();
    while($temp[7] =~ /(cpd\d{5})\[(\w)\]/g){
	$Cmpts{$2}{$1}=1;
    }
    next unless scalar(keys %Cmpts)>1;

    $Main_Cmpts{$temp[0]}=\%Cmpts;
}
close(FH);

open(FH, "< ../Biochemistry/reactions.master.tsv");
open(OUT, "> Transport_Main_Compartment_Translation.tsv");
$header=1;
while(<FH>){
    chomp;
    if($header){$header--;next}
    @temp=split(/\t/,$_,-1);

    my %Cmpts=();
    while($temp[6] =~ /(cpd\d{5})\[(\w)\]/g){
	$Cmpts{$1}=$2;
    }
    next unless exists($Main_Cmpts{$temp[0]});

    my $GenCmpt=-1;
    my $OrgCmpt=-1;
    if(exists($Main_Cmpts{$temp[0]}{'c'})){
	foreach my $cpd (keys %{$Main_Cmpts{$temp[0]}{'c'}}){
	    $GenCmpt=$Cmpts{$cpd} if exists($Cmpts{$cpd});
	}
	$GenCmpt = 0 if $GenCmpt == -1 && scalar(keys %Cmpts)==1;
	$OrgCmpt='c';
    }elsif(exists($Main_Cmpts{$temp[0]}{'p'})){
	foreach my $cpd (keys %{$Main_Cmpts{$temp[0]}{'p'}}){
	    $GenCmpt=$Cmpts{$cpd} if exists($Cmpts{$cpd});
	}
	$GenCmpt = 0 if $GenCmpt == -1 && scalar(keys %Cmpts)==1;
	$OrgCmpt='p';
    }elsif(exists($Main_Cmpts{$temp[0]}{'d'})){
	foreach my $cpd (keys %{$Main_Cmpts{$temp[0]}{'d'}}){
	    $GenCmpt=$Cmpts{$cpd} if exists($Cmpts{$cpd});
	}
	$GenCmpt = 0 if $GenCmpt == -1 && scalar(keys %Cmpts)==1;
	$OrgCmpt='d';
    }

    if($GenCmpt eq "-1"){
	print scalar(keys %Cmpts),"\n";
    }

    print OUT $temp[0],"\t",$OrgCmpt,"\t",$GenCmpt,"\n";
}
close(FH);
close(OUT);

__END__

#The script couldn't handle a small number of reactions due to differing compound ids
#These were manually edited in the output file:
#rxn13439        c       1
#rxn13645        c       0
#rxn13936        c       0
