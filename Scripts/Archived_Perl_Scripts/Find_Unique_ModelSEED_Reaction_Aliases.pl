#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my $header=1;
my @temp=();
my $JSON=undef;

my %Aliases_MS=();
my %MS_Aliases=();
foreach my $type ("Compounds","Reactions"){
    open(FH, "< Provenance/${type}_Aliases.tsv");
    $header=1;
    while(<FH>){
	chomp;
	if($header){$header--;next}
	@temp=split(/\t/,$_);

	if($temp[0] =~ /(cpd|rxn)/){
	    foreach my $id (split(/\|/,$temp[0])){
		$Aliases_MS{$temp[2]}{$id}=1;
		$MS_Aliases{$id}{$temp[2]}{$temp[3]}=1;
	    }
	}elsif($temp[1] =~ /(cpd|rxn)/){
	    foreach my $id (split(/\|/,$temp[1])){
		$Aliases_MS{$temp[2]}{$id}=1;
		$MS_Aliases{$id}{$temp[2]}{$temp[3]}=1;
	    }
	}else{
	    print $_,"\n";
	}
    }
    close(FH);
}

open(FH, "< $ENV{SEAVER_PROJECT}ModelSEEDDatabase/Biochemistry/reactions.json");
$JSON="";while(<FH>){$JSON.=$_;}
close(FH);
my %Biochemical_Reactions=%{from_json($JSON)};

my %Unique_MS_Aliases=();
my %Unique_Aliases_MS=();
foreach my $rxn (sort keys %Biochemical_Reactions){
    foreach my $alias (keys %{$MS_Aliases{$rxn}}){
	$Unique_MS_Aliases{$rxn}{$alias}={};
	$Unique_Aliases_MS{$alias}{$rxn}={};
    }

    foreach my $link (split(/;/,$Biochemical_Reactions{$rxn}{linked_reaction})){
	foreach my $alias (keys %{$MS_Aliases{$rxn}}){
	    $Unique_MS_Aliases{$link}{$alias}={};
	    $Unique_Aliases_MS{$alias}{$link}={};
	}

	foreach my $alias (keys %{$MS_Aliases{$link}}){
	    $Unique_MS_Aliases{$rxn}{$alias}={};
	    $Unique_Aliases_MS{$alias}{$rxn}={};
	    $Unique_MS_Aliases{$link}{$alias}={};
	    $Unique_Aliases_MS{$alias}{$link}={};
	}
    }
}

foreach my $rxn (sort keys %Biochemical_Reactions){
    foreach my $alias (keys %{$MS_Aliases{$rxn}}){
	foreach my $source (keys %{$MS_Aliases{$rxn}{$alias}}){
	    $Unique_MS_Aliases{$rxn}{$alias}{$source}=1;
	}
    }

    foreach my $link (split(/;/,$Biochemical_Reactions{$rxn}{linked_reaction})){
	foreach my $alias (keys %{$MS_Aliases{$rxn}}){
	    foreach my $source (keys %{$MS_Aliases{$rxn}{$alias}}){
		$Unique_MS_Aliases{$link}{$alias}{$source}=1;
	    }
	}

	foreach my $alias (keys %{$MS_Aliases{$link}}){
	    foreach my $source (keys %{$MS_Aliases{$link}{$alias}}){
		$Unique_MS_Aliases{$link}{$alias}{$source}=1;
		$Unique_MS_Aliases{$rxn}{$alias}{$source}=1;
	    }
	}
    }
}

open(OUT, "> Unique_ModelSEED_Reaction_Aliases.txt");
print OUT "ModelSEED ID\tExternal ID\tSource\n";
foreach my $rxn (sort keys %Biochemical_Reactions){
    foreach my $alias (sort keys %{$Unique_MS_Aliases{$rxn}}){
	print OUT $rxn,"\t",$alias,"\t",join("|",sort(keys %{$Unique_MS_Aliases{$rxn}{$alias}})),"\n";
    }
}
close(OUT);
