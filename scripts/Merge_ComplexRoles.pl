#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header=1;

######################
#Read in current roles
######################
open(FH, "< ../Mappings/probmodelseed-mapping/Mapping_Roles.txt");
my @headers = split(/\t/,<FH>);
chomp($headers[$#headers]);
my %Current_Role_Names = ();
my %Current_Role_IDs = ();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $Current_Role_Names{$temp[1]}=$temp[0];
    $Current_Role_IDs{$temp[0]}=$temp[1];
}
close(FH);

open(FH, "< ../Mappings/ComplexRoles.kegg.tsv");
my %Transform_IDs=();
my %New_Role_IDs = ();
$header = 1;
while(<FH>){
    chomp;
    if($header){$header--;next}
    @temp=split(/\t/,$_,-1);

    my $Is_Present=0;
    if(exists($Current_Role_Names{$temp[5]})){
	$Transform_IDs{$temp[4]}=$Current_Role_Names{$temp[5]};
	$Is_Present=1;
    }

    #attempt aliases just in case
    foreach my $set_alias (split(/;/,$temp[8])){
	my ($set,$alias)=split(/:/,$set_alias);
	next if !$alias;

	if(exists($Current_Role_Names{$alias})){

	    #WARNING: Conflicting names and ids exist
	    if(exists($Transform_IDs{$temp[4]}) && $Transform_IDs{$temp[4]} ne $Current_Role_Names{$alias}){
		print "Warning: Conflicting ids/names\n";
		next;
	    }

	    $Transform_IDs{$temp[4]}=$Current_Role_Names{$alias};
	    $Is_Present=1;

	}
    }

    next if $Is_Present;
    $New_Role_IDs{$temp[4]}=$temp[5];
}
close(FH);

open(OUT, "> ../Mappings/probmodelseed-mapping/Mapping_Roles.txt");
print OUT join("\t",@headers),"\n";
foreach my $id (sort keys %Current_Role_IDs){
    print OUT $id,"\t",$Current_Role_IDs{$id},"\tnull\n";
}
foreach my $id (sort keys %New_Role_IDs){
    print OUT $id,"\t",$New_Role_IDs{$id},"\tnull\n";
}
close(OUT);
