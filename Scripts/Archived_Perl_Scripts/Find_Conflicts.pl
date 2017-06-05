#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my @temp=();
my %Source_Conflicts=();

my $BPath = "../../Biochemistry/";
open(FH, "< ".$BPath."Aliases/Compounds_Aliases.tsv");
while(<FH>){
    chomp;
    my ($pmsid,$dmsid,$sid,$source)=split(/\t/,$_,-1);
    next unless $source =~ /^(KEGG|BiGG)$/;

    $Source_Conflicts{$source}{"compounds"}{$sid}{'p'}{$pmsid}=1;
    $Source_Conflicts{$source}{"compounds"}{$sid}{'d'}{$dmsid}=1;
}
close(FH);

open(FH, "< ".$BPath."Aliases/Reactions_Aliases.tsv");
while(<FH>){
    chomp;
    my ($pmsid,$dmsid,$sid,$source)=split(/\t/,$_,-1);
    next unless $source =~ /^(KEGG|BiGG)$/;

    $Source_Conflicts{$source}{"reactions"}{$sid}{'p'}{$pmsid}=1;
    $Source_Conflicts{$source}{"reactions"}{$sid}{'d'}{$dmsid}=1;
}
close(FH);

foreach my $source (sort keys %Source_Conflicts){
    foreach my $type (sort keys %{$Source_Conflicts{$source}}){
	open(OUT, "> ".$source."_".$type."_Conflicts.tsv");
	print OUT "Source\tType\tExternal ID\tNew ModelSEED ID\tOld ModelSEED ID\n";
	foreach my $id (sort keys %{$Source_Conflicts{$source}{$type}}){
	    my $pstring = join("|", sort keys %{$Source_Conflicts{$source}{$type}{$id}{'p'}});
	    my $dstring = join("|", sort keys %{$Source_Conflicts{$source}{$type}{$id}{'d'}});
	    if($dstring && $pstring ne $dstring){
		print OUT $source,"\t",$type,"\t",$id,"\t",$pstring,"\t",$dstring,"\n";
	    }
	}
	close(OUT);
    }
}
