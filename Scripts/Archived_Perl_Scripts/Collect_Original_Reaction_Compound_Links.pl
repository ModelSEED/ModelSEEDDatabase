#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();

#Remember that KEGG and MetaCyc aliases have artificial ones inserted by you
opendir(my $dh, "../");
my @Files = grep { $_ =~ /_Aliases\.txt$/ } readdir($dh);
closedir($dh);

my %Sources=();
foreach my $file (@Files){
    open(FH, "< ../${file}");
    while(<FH>){
	chomp;
	@temp=split(/\t/,$_,-1);
	$Sources{$temp[2]}{$temp[1]}=1;
    }
}

opendir($dh, $ENV{SEAVER_PROJECT}."Model_Import_Files/");
my @Models = grep { $_ =~ /-reactions\.tbl/ && $_ !~ /Cmpt/ } readdir($dh);
closedir($dh);

open(OUT, "> Original_Reaction_Compound_Links.tsv");
foreach my $Model (@Models){
    my %Reaction_Links=();
    my $model = (split(/-/,$Model))[0];
    $model =~ s/import//;

    if($model =~ /(KEGG).+/){
	$model = $1;
    }

    if(exists($Sources{$model})){
	open(FH, "< ".$ENV{SEAVER_PROJECT}."Model_Import_Files/".$Model);
	my $header = 1;
	my $equation_header=0;
	while(<FH>){
	    chomp;
	    @temp=split(/\t/,$_,-1);
	    if(scalar(@temp)<3){
		@temp=split(/;/,$_,-1);
	    }
	    if($header){
		for(my $i=0;$i<scalar(@temp);$i++){
		    if($temp[$i] eq "EQUATION"){
			$equation_header=$i;
		    }
		}
		$header--;
		next;
	    }

	    #parse equation
	    my @equation = split(/\s+/,$temp[$equation_header]);
	    foreach my $entity (@equation){
		$entity =~ s/\[\w\]$//;
		if(exists($Sources{$model}{$entity})){
		    $Reaction_Links{$temp[0]}{$entity}=1;
		}
	    }
	}
	close(FH);

	foreach my $rxn (sort keys %Reaction_Links){
	    print OUT $model,"\t",$rxn,"\t",join("|",sort keys %{$Reaction_Links{$rxn}}),"\n";
	}
    }
}

foreach my $model ("KEGG","MetaCyc"){
    my %Reaction_Links=();
    open(FH, "< Primary_Databases/".$model."_Reactions.tbl");
    my $header = 1;
    my $equation_header=0;
    while(<FH>){
	chomp;
	@temp=split(/\t/,$_,-1);
	if(scalar(@temp)<3){
	    @temp=split(/;/,$_,-1);
	}
	if($header){
	    for(my $i=0;$i<scalar(@temp);$i++){
		if($temp[$i] eq "EQUATION"){
		    $equation_header=$i;
		}
	    }
	    $header--;
	    next;
	}
	
	#parse equation
	my @equation = split(/\s+/,$temp[$equation_header]);
	foreach my $entity (@equation){
	    $entity =~ s/\[\w\]$//;
	    if(exists($Sources{$model}{$entity})){
		$Reaction_Links{$temp[0]}{$entity}=1;
	    }
	}
    }
    close(FH);

    foreach my $rxn (sort keys %Reaction_Links){
	print OUT $model,"\t",$rxn,"\t",join("|",sort keys %{$Reaction_Links{$rxn}}),"\n";
    }
}
close(OUT);
