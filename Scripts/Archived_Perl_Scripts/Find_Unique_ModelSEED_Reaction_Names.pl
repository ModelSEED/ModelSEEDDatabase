#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my $header=1;
my @temp=();
my $JSON=undef;

my %Aliases_MS=();
my %MS_Aliases=();
my %Searchnames=();
foreach my $type ("Compounds","Reactions"){
    open(FH, "< Provenance/Names_${type}_Aliases.tsv");
    $header=1;
    while(<FH>){
	chomp;
	if($header){$header--;next}
	@temp=split(/\t/,$_);

	if($temp[3] eq 'searchname'){
	    $Searchnames{$temp[2]}=1;
	    next;
	}

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

open(FH, "< Source_Classifiers.txt");
my %Models = ();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    if($temp[1] =~ /Model/){
	$Models{$temp[0]}=1;
    }
}
close(FH);

open(OUT, "> Unique_ModelSEED_Reaction_Names.txt");
print OUT "ModelSEED ID\tExternal ID\tSource\n";
foreach my $rxn (sort keys %Biochemical_Reactions){
    foreach my $alias (sort keys %{$Unique_MS_Aliases{$rxn}}){
	next if $alias eq "-";
	next if length($alias)<2;

	my $uc_alias = uc(substr($alias,0,1)).substr($alias,1);
	next if exists($Searchnames{$alias}) && exists($Unique_MS_Aliases{$rxn}{$uc_alias});

	my $is_model=0;
	foreach my $source (keys %{$MS_Aliases{$rxn}{$alias}}){
	    $Unique_MS_Aliases{$rxn}{$alias}{$source}=1;
	    if(exists($Models{$source})){
		$is_model=1;
	    }
	}

	my $mod_alias=$alias;
	if($is_model){
	    #remove metabolite 'M_' prefix
	    $mod_alias =~ s/^M_//;

	}
#	print $mod_alias,"\n" if $is_model;

	print OUT $rxn,"\t",$alias,"\t",join("|",sort(keys %{$Unique_MS_Aliases{$rxn}{$alias}})),"\n";
    }
}
close(OUT);
