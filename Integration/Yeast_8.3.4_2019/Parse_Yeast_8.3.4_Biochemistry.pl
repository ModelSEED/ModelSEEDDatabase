#!/usr/bin/env perl
#use PlantSEED::Formulas;
use warnings;
use strict;
my @temp=();

my $Compounds = "yeastGEM_Compound_Table.txt";
open(FH, "< $Compounds");
my $header=1;
my %Original_Compounds=();
while(<FH>){
    chomp;
    if($header){$header--;next}
    @temp=split(/\t/,$_);

    my $cpd_cpt = $temp[0];

    #Convert ASCII codes
    $temp[0] =~ s/__91__/[/;
    $temp[0] =~ s/__93__/]/;

    #Clean up identifier
    $temp[0] =~ s/^s_+//;

    #Remove Compartment
    my $cpd = $temp[0];
    $cpd =~ s/\[\w+\]$//;

    #Strip name
    $temp[1] =~ s/(\s\[[\w\s]+\])$//;

    #Define KEGG
    $temp[5] = "" if !$temp[5];

    $Original_Compounds{$cpd_cpt}={'ID'=>$cpd,
				   'NAMES'=>$temp[1],
				   'KEGG'=>$temp[5]};
}
close(FH);

my $filestub = $Compounds;
$filestub =~ s/_Compound_Table\.txt$//;

open(OUT, "> ".$filestub."_Compounds.tbl");
my @Headers=("ID","NAMES","KEGG");
print OUT join("\t",@Headers),"\n";
foreach my $id (sort keys %Original_Compounds){
    foreach my $h (@Headers){
	print OUT $Original_Compounds{$id}{$h};
	print OUT "\t" unless $h eq $Headers[$#Headers];
    }
    print OUT "\n";
}
close(OUT);

__END__

open(OUT, "> ".$Output_Root."KEGG_Glycans.tbl");
print OUT "ID\t".join("\t",@Headers)."\n";
foreach my $id (sort keys %$Glycans){
    print OUT $id,"\t";
    foreach my $h (@Headers){
	print OUT $Glycans->{$id}{$h};
	print OUT "\t" unless $h eq $Headers[$#Headers];
    }
    print OUT "\n";
}
close(OUT);

open(OUT, "> ".$Output_Root."KEGG_Reactions.tbl");
@Headers = ("NAMES","KEGG","EQUATION","ECs");
print OUT "ID\t".join("\t",@Headers)."\n";
foreach my $c (sort { $a <=> $b } keys %$Reactions){
    my @ids=sort grep { substr($_,0,1) eq "R" } split /\|/,$Reactions->{$c}{"KEGG"};

    print OUT $ids[0],"\t";
    foreach my $h (@Headers){
	if(!defined($Reactions->{$c}{$h})){
	    print STDERR "Missing $h\n";
	}
	print OUT $Reactions->{$c}{$h};
	print OUT "\t" unless $h eq $Headers[$#Headers];
    }
    print OUT "\n";
}
close(OUT);
