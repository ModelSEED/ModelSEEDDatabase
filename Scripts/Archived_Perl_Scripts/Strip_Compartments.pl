#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();

#For searchname
use Bio::KBase::ObjectAPI::KBaseBiochem::Compound;

open(FH, "< Names_Reactions_Aliases.tsv");
my %Lines=();
my $Header=<FH>;chomp $Header;
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next unless $temp[-1] eq "name";
    $temp[-2] =~ s/periplasm-cytoplasm //;
    $temp[-2] =~ s/\(extracellular to periplasm\)//;
    $temp[-2] =~ s/\(cytoplasm to periplasm\)//;
    $temp[-2] =~ s/, cytoplasm to periplasm//;
    $temp[-2] =~ s/\[periplasm-cytosol\]//;
    $temp[-2] =~ s/\[extraorganism-periplasm\]//;
    $temp[-2] =~ s/\(periplasm to extracellular\)//;
    $temp[-2] =~ s/\(periplasm to cytoplasm\)//;
    $temp[-2] =~ s/periplasm to cytoplasm, //;

    $temp[-2] =~ s/\(periplasm\)//;
    $temp[-2] =~ s/\(periplam\)//;
    $temp[-2] =~ s/\(to periplasm\)//;
    $temp[-2] =~ s/\( periplasm \)//;
    $temp[-2] =~ s/\( periplasm\)//;
    $temp[-2] =~ s/\(periplasmic\)//;
    $temp[-2] =~ s/periplasmic, //;
    $temp[-2] =~ s/, periplasm//;
    $temp[-2] =~ s/periplasm, //;
    $temp[-2] =~ s/periplasmic//;
    $temp[-2] =~ s/ periplasm//;
    $temp[-2] =~ s/periplasm//;

    $temp[-2] =~ s/\(Plastidial membrane\)//;
    $temp[-2] =~ s/\(Mitochondrial\/Plastidial membrane\)//;
    $temp[-2] =~ s/\(Plastidal membrane\)//;
    $temp[-2] =~ s/\(Plastid membrane\)//;
    $temp[-2] =~ s/\(plastidal\)//;

    $temp[-2] =~ s/ \(chloroplast\)//;
    $temp[-2] =~ s/ \(into chloroplast\)//;
    $temp[-2] =~ s/, chloroplast$//;
    $temp[-2] =~ s/, chloroplastic$//;
    $temp[-2] =~ s/_between_cytosol_and_plastid//;
    $temp[-2] =~ s/ chloroplast \///;
    $temp[-2] =~ s/ to chloroplast //;
    $temp[-2] =~ s/, chloroplast //;
    $temp[-2] =~ s/ chloroplast,//;
    $temp[-2] =~ s/chloroplast transport/transport/;
    $temp[-2] =~ s/^anaplastic//;
    $temp[-2] =~ s/^Chlorplastic//;
    $temp[-2] =~ s/^Chloroplastice//;
    $temp[-2] =~ s/^Chlroroplastic//;
    $temp[-2] =~ s/^Chloroplast//i;

    $temp[-2] =~ s/mitochondrial transport/transport/;
    $temp[-2] =~ s/mtiochondrial transport/transport/;
    $temp[-2] =~ s/mitochondrial  transport/transport/;
    $temp[-2] =~ s/\(mitochondria\)//;
    $temp[-2] =~ s/\(mitochondrial\)//i;
    $temp[-2] =~ s/\(Mitochondrial membrane\)//;
    $temp[-2] =~ s/mitochondrial electron transport/electron transport/;
    $temp[-2] =~ s/out of mitochondria //;
    $temp[-2] =~ s/ to mitochondriat//;
    $temp[-2] =~ s/ from mitochondia to cytoplasm//;
    $temp[-2] =~ s/, cytosolic\/mitochondrial//;
    $temp[-2] =~ s/ from mitochondria to cytosol//;
    $temp[-2] =~ s/_between_cytosol_and_mitochondrion//;
    $temp[-2] =~ s/ into mitochondria//;
    $temp[-2] =~ s/ to mitochondria//;
    $temp[-2] =~ s/, mitochondrial 1//;
    $temp[-2] =~ s/-mitochondrial//;
    $temp[-2] =~ s/, mitochondrial,/,/;
    $temp[-2] =~ s/, mitochondria$//;
    $temp[-2] =~ s/, mitochondrial$//;
    $temp[-2] =~ s/, mitochondiral$//;
    $temp[-2] =~ s/, mitochodnrial$//;
    $temp[-2] =~ s/ mitochondrial / /;
    $temp[-2] =~ s/\(mitochondrial / /;
    $temp[-2] =~ s/^Mitochindrial//;
    $temp[-2] =~ s/^Mitocondrial//;
    $temp[-2] =~ s/^Mitochondrial//i;

    $temp[-2] =~ s/, cytoplasmic,/,/;
    $temp[-2] =~ s/, cytoplasmic/,/g;
    $temp[-2] =~ s/ cytoplasmic$//g;
    $temp[-2] =~ s/ exocytoplasmic//g;
    $temp[-2] =~ s/\[extraorganism-cytosol\]//;
    $temp[-2] =~ s/\[cytosol-extraorganism\]//;
    $temp[-2] =~ s/_between_cytosol_and_peroxisome//;
    $temp[-2] =~ s/\(cytosol to extracellular\)//;
    $temp[-2] =~ s/\(cytosol to nucleus\)//;
    $temp[-2] =~ s/, glycosome to cytosol//;
    $temp[-2] =~ s/, cytosol to extracellular//;
    $temp[-2] =~ s/\(extracellular to cytosol\)//;
    $temp[-2] =~ s/, located in cytosol//;
    $temp[-2] =~ s/, cytosolic$//;
    $temp[-2] =~ s/, cytosolic / /;
    $temp[-2] =~ s/, cytosol$//;
    $temp[-2] =~ s/^cytosolic //i;
    $temp[-2] =~ s/^cytosol //i;

    $temp[-2] =~ s/\(extracellular\)//;
    $temp[-2] =~ s/\(extracellular$//;
    $temp[-2] =~ s/, extracellular$//;
    $temp[-2] =~ s/ extracellular / /;
    $temp[-2] =~ s/\(extracellular membrane\)//;
    $temp[-2] =~ s/\(to extracellular\)//;
    $temp[-2] =~ s/^extracellular//i;
    $temp[-2] =~ s/^extracellualr//i;

    $temp[-2] =~ s/, endoplasmic reticular$//;
    $temp[-2] =~ s/, endoplamic reticular$//;
    $temp[-2] =~ s/, endoplasmic reticular / /;
    $temp[-2] =~ s/, endoplasmic reticulum$//;
    $temp[-2] =~ s/endoplasmic reticular transport/transport/;
    $temp[-2] =~ s/endoplasmic reticulum transport/transport/;
    $temp[-2] =~ s/endoplamic reticular transport/transport/;
    $temp[-2] =~ s/endoplamic reticulum transport/transport/;
    $temp[-2] =~ s/sarco\(endo\)plasmic reticulum//;
    $temp[-2] =~ s/sarcoplasmic reticulum//;

    $temp[-2] =~ s/, golgi$//i;
    $temp[-2] =~ s/, golgi apparatus$//i;
    $temp[-2] =~ s/\(Golgi\)//;
    $temp[-2] =~ s/\(Golgi apparatus\)//;
    $temp[-2] =~ s/Golgi transport/transport/i;
    $temp[-2] =~ s/^Golgi//i;

    $temp[-2] =~ s/, nuclear$//;
    $temp[-2] =~ s/, nucleus$//;
    $temp[-2] =~ s/nuclear transport/transport/;
    $temp[-2] =~ s/, nuclear//;
    $temp[-2] =~ s/ nuclear / /;
    $temp[-2] =~ s/^nuclear//i;

    $temp[-2] =~ s/glyoxysomal transport/transport/;
    $temp[-2] =~ s/, glyoxysomal$//;
    $temp[-2] =~ s/, glyoxysome$//;
    $temp[-2] =~ s/, glyxoxysome$//;
    $temp[-2] =~ s/, glyxosyome$//;

    $temp[-2] =~ s/\(Plasma membrane\)//;
    $temp[-2] =~ s/, plasma membrane$//;
    $temp[-2] =~ s/ plasma membrane / /;
    $temp[-2] =~ s/ outer membrane / /i;
    $temp[-2] =~ s/\(outer membrane /(/i;
    $temp[-2] =~ s/ outher membrane / /i;
    $temp[-2] =~ s/^plasma membrane / /i;
    $temp[-2] =~ s/\(Peroxisomal membrane\)//;
    $temp[-2] =~ s/\(Vacuolar membrane\)//;
    $temp[-2] =~ s/^membrane[\s-]//i;

    $temp[-2] =~ s/peroxisomal transport/transport/;
    $temp[-2] =~ s/peroxisomal shuttle/shuttle/;
    $temp[-2] =~ s/, peroxisomal$//;
    $temp[-2] =~ s/, peroxisome$//;
    $temp[-2] =~ s/, peroxisomal,/,/;
    $temp[-2] =~ s/into peroxisome$//;
    $temp[-2] =~ s/out of peroxisome$//;
    $temp[-2] =~ s/\(peroxisomal\)//i;
    $temp[-2] =~ s/^peroxisomal//i;
    $temp[-2] =~ s/^peroxisome//i;

    $temp[-2] =~ s/, vacuolar$//i;
    $temp[-2] =~ s/, vacuole$//i;
    $temp[-2] =~ s/\(vacuolar\)//i;
    $temp[-2] =~ s/vacuolar transport/transport/;
    $temp[-2] =~ s/vacuolar \'transport\'/transport/;
    $temp[-2] =~ s/^vacuolar//i;

    $temp[-2] =~ s/__/_/;
    $temp[-2] =~ s/  / /;
    $temp[-2] =~ s/^\s+//;
    $temp[-2] =~ s/\s+$//;

    foreach my $new_rxn (split(/\|/,$temp[0])){
	$Lines{$temp[-2]}{'new'}{$new_rxn}=1;
    }

    foreach my $old_rxn (split(/\|/,$temp[1])){
	$Lines{$temp[-2]}{'old'}{$old_rxn}=1;
    }
}
close(FH);

open(OUT, "> Names_Reactions_Aliases.tmp");
print OUT $Header."\n";
foreach my $name (sort keys %Lines){
    print OUT join("|",sort keys %{$Lines{$name}{'new'}}),"\t";
    print OUT join("|",sort keys %{$Lines{$name}{'old'}}),"\t";
    print OUT $name,"\tname\n";
}

foreach my $name (sort keys %Lines){
    print OUT join("|",sort keys %{$Lines{$name}{'new'}}),"\t";
    print OUT join("|",sort keys %{$Lines{$name}{'old'}}),"\t";
    print OUT Bio::KBase::ObjectAPI::KBaseBiochem::Compound::nameToSearchname($name),"\tsearchname\n";
}

close(OUT);
