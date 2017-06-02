#!/usr/bin/env perl
use warnings;
use strict;
use Pathways;
my @temp=();

if(!$ARGV[0]){
    print STDERR "Usage: ./Create_Pathway_ReactionSets.pl <KEGG|MetaCyc|AraCyc|etc.>\n";
    print join("\n", Pathways::list_cyc_pathways()),"\n";
    exit;
}

my $DB=$ARGV[0];
my $Pathways=Pathways::get_pathways({db=>$DB});

open(OUT, "> ../Pathways/".$DB.".pathways");
print OUT "Source\tSource ID\tName\tAliases\tReactions\n";
foreach my $Pathway (sort keys %$Pathways){
    print OUT $DB,"\t",$Pathway,"\t",$Pathways->{$Pathway}{"Name"},"\t",$Pathways->{$Pathway}{"Class"},"\t",join("|", @{$Pathways->{$Pathway}{"Reactions"}}),"\n";
}
close(OUT);
