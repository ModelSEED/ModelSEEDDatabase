#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my @temp=();
my $header=1;
my %Source_Conflicts=();

use lib '../../Libs/Perl/';
use InChIs;
use Formulas;

my $BPath = "../../Biochemistry/";
open(FH, "< ".$BPath."Structures/Structures.tsv");
my %Source_InChIs=();
my %Source_Formulas=();
my %Source_Atoms=();
my %Source_Layers=();
while(<FH>){
    chomp;
    my ($msid,$id,$source,$inchi)=split(/\t/,$_,-1);
    if($header){$header--;next}
    my ($formula,$layers) = InChIs::parse($inchi);
    my $new_inchi = InChIs::build($formula,$layers,{remove=>{'p'=>1}});
    $Source_InChIs{$msid}{$new_inchi}{$id}=1;
    $Source_Formulas{$msid}{$formula}{$id}=1;
    my $atoms = Formulas::parse($formula);
    foreach my $atom ( sort keys %$atoms ){
	$Source_Atoms{$msid}{$atom}{$atoms->{$atom}}{$id}=1;
    }
    foreach my $layer ( grep { $_ ne 'p' } sort keys %$layers){
	next if !$layers->{$layer};
	$Source_Layers{$msid}{$layer}{$layers->{$layer}}{$id}=1;
    }
}
close(FH);

open(OUT, "> InChI_Formula_Conflicts.tsv");
print OUT "ModelSEED ID\tFormula\tExternal IDs\n";
my %Skip_Formulas=();
foreach my $msid (sort grep { scalar(keys %{$Source_Formulas{$_}})>1 } keys %Source_Formulas){
    $Skip_Formulas{$msid}=1;
    foreach my $formula (sort keys %{$Source_Formulas{$msid}}){
	print OUT $msid,"\t",$formula,"\t",join("|",sort keys %{$Source_Formulas{$msid}{$formula}}),"\n";
    }
}
close(OUT);

open(OUT, "> InChI_Atoms_Conflicts.tsv");
print OUT "ModelSEED ID\tAtom\tExternal IDs\n";
foreach my $msid ( grep { exists($Skip_Formulas{$_}) } sort keys %Source_Atoms){
    foreach my $atom (sort grep { scalar(keys %{$Source_Atoms{$msid}{$_}})>1 } keys %{$Source_Atoms{$msid}}){
	foreach my $stoichiometry (sort keys %{$Source_Atoms{$msid}{$atom}}){
	    print OUT $msid,"\t",$atom,":",$stoichiometry,"\t",join("|",sort keys %{$Source_Atoms{$msid}{$atom}{$stoichiometry}}),"\n";
	}
    }
}
close(OUT);

open(OUT, "> InChI_Layers_Conflicts.tsv");
print OUT "ModelSEED ID\tLayer\tExternal IDs\n";
foreach my $msid ( grep { !exists($Skip_Formulas{$_}) } sort keys %Source_Layers){
    foreach my $layer (sort grep { scalar(keys %{$Source_Layers{$msid}{$_}})>1 } keys %{$Source_Layers{$msid}}){
	foreach my $layer_contents (sort keys %{$Source_Layers{$msid}{$layer}}){
	    print OUT $msid,"\t",$layer,":",$layer_contents,"\t",join("|",sort keys %{$Source_Layers{$msid}{$layer}{$layer_contents}}),"\n";
	}
    }
}
close(OUT);

open(OUT, "> InChI_Full_InChI_Conflicts.tsv");
print OUT "ModelSEED ID\tInChI\tExternal IDs\n";
foreach my $msid (sort grep { scalar(keys %{$Source_InChIs{$_}})>1 } keys %Source_InChIs){
    foreach my $inchi (sort keys %{$Source_InChIs{$msid}}){
	print OUT $msid,"\t",$inchi,"\t",join("|",sort keys %{$Source_InChIs{$msid}{$inchi}}),"\n";
    }
}
close(OUT);
