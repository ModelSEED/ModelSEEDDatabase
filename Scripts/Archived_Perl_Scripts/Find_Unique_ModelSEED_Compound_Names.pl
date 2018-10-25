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
    open(FH, "< Provenance/Names_${type}_Aliases.tsv");
    $header=1;
    while(<FH>){
	chomp;
	if($header){$header--;next}
	@temp=split(/\t/,$_);

	next if $temp[3] eq "searchname";

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

open(FH, "< $ENV{SEAVER_PROJECT}ModelSEEDDatabase/Biochemistry/compounds.json");
$JSON="";while(<FH>){$JSON.=$_;}
close(FH);
my %Biochemical_Compounds=%{from_json($JSON)};

my %Unique_MS_Aliases=();
my %Unique_Aliases_MS=();
foreach my $cpd (sort keys %Biochemical_Compounds){
    foreach my $alias (keys %{$MS_Aliases{$cpd}}){
	$Unique_MS_Aliases{$cpd}{$alias}={};
	$Unique_Aliases_MS{$alias}{$cpd}={};
    }

    foreach my $link (split(/;/,$Biochemical_Compounds{$cpd}{linked_compound})){
	foreach my $alias (keys %{$MS_Aliases{$cpd}}){
	    $Unique_MS_Aliases{$link}{$alias}={};
	    $Unique_Aliases_MS{$alias}{$link}={};
	}

	foreach my $alias (keys %{$MS_Aliases{$link}}){
	    $Unique_MS_Aliases{$cpd}{$alias}={};
	    $Unique_Aliases_MS{$alias}{$cpd}={};
	    $Unique_MS_Aliases{$link}{$alias}={};
	    $Unique_Aliases_MS{$alias}{$link}={};
	}
    }
}

foreach my $cpd (sort keys %Biochemical_Compounds){
    foreach my $alias (keys %{$MS_Aliases{$cpd}}){
	foreach my $source (keys %{$MS_Aliases{$cpd}{$alias}}){
	    $Unique_MS_Aliases{$cpd}{$alias}{$source}=1;
	}
    }

    foreach my $link (split(/;/,$Biochemical_Compounds{$cpd}{linked_compound})){
	foreach my $alias (keys %{$MS_Aliases{$cpd}}){
	    foreach my $source (keys %{$MS_Aliases{$cpd}{$alias}}){
		$Unique_MS_Aliases{$link}{$alias}{$source}=1;
	    }
	}

	foreach my $alias (keys %{$MS_Aliases{$link}}){
	    foreach my $source (keys %{$MS_Aliases{$link}{$alias}}){
		$Unique_MS_Aliases{$link}{$alias}{$source}=1;
		$Unique_MS_Aliases{$cpd}{$alias}{$source}=1;
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

open(OUT, "> Unique_ModelSEED_Compound_Names.txt");
print OUT "ModelSEED ID\tExternal ID\tSource\n";
foreach my $cpd (sort keys %Biochemical_Compounds){
    foreach my $alias (sort keys %{$Unique_MS_Aliases{$cpd}}){

	my $is_model=0;
	foreach my $source (keys %{$MS_Aliases{$cpd}{$alias}}){
	    $Unique_MS_Aliases{$cpd}{$alias}{$source}=1;
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

	print OUT $cpd,"\t",$alias,"\t",join("|",sort(keys %{$Unique_MS_Aliases{$cpd}{$alias}})),"\n";
    }
}
close(OUT);
