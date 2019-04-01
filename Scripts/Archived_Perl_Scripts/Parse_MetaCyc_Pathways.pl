#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();

my $Database_Root = "/homes/seaver/Biochemistry_Mirrors/bioinformatics.ai.sri.com/metacyc/dist/22.5/data/";
my $Structure_Root = "/homes/seaver/Projects/ModelSEEDDatabase/Biochemistry/Structures/MetaCyc/";
my $Output_Root = "/homes/seaver/Projects/ModelSEEDDatabase/Biochemistry/Aliases/Provenance/Primary_Databases/";

my %PIDsNames=();
my %PIDsTypes=();
my %IsSuper=();
open(FH, "< ".$Database_Root."/classes.dat");
my ($Field,$Data)=("","");
my $ID='';
while(<FH>){
    chomp;
    if(substr($_,0,1) eq "#"){next;}
    @temp=split/\s-\s/;
    
    $Field=$temp[0];
    $Data=join("-",@temp[1..$#temp]);
    next unless $Field;

    next if $Data eq "FRAMES" || $Data eq "Generalized-Reactions" || $Data eq "Pathways" || $Data eq "Super-Pathways";
    
    $Field=~s/\s*$//;
    $Data=~s/^\s*//;
    $Data=~s/\s*$//;
    
    if($Field eq "UNIQUE-ID"){
	$ID=$Data;
    }elsif($Field eq "COMMON-NAME"){
	$PIDsNames{$ID}{$Data}=1;
    }elsif($Field eq "TYPES"){
	$PIDsTypes{$ID}{$Data}=1;
	$IsSuper{$Data}=1;
    }
}
close(FH);

my %ReactionList=();
my %PIDsSuper=();
open(FH, "< ".$Database_Root."/pathways.dat");
($Field,$Data)=("","");
$ID='';
while(<FH>){
    chomp;
    if(substr($_,0,1) eq "#"){next;}
    @temp=split/\s-\s/;
    
    $Field=$temp[0];
    $Data=join("-",@temp[1..$#temp]);
    next unless $Field;
    
    $Field=~s/\s*$//;
    $Data=~s/^\s*//;
    $Data=~s/\s*$//;
    
    next if $Data eq "FRAMES" || $Data eq "Generalized-Reactions" || $Data eq "Pathways" || $Data eq "Super-Pathways";

    if($Field eq "UNIQUE-ID"){
	#07/25
	#All UNIQUE-IDs are UNIQUE
	$ID=$Data;
    }elsif($Field eq "COMMON-NAME"){
	$PIDsNames{$ID}{$Data}=1;
    }elsif($Field eq "TYPES"){
	$PIDsTypes{$ID}{$Data}=1;
	$IsSuper{$Data}=1;
    }elsif($Field eq "REACTION-LIST"){
	$ReactionList{$ID}{$Data}=1;
    }elsif($Field eq "SUPER-PATHWAYS"){
	$PIDsSuper{$ID}{$Data}=1;
	$IsSuper{$Data}=1;
    }elsif($Field eq "IN-PATHWAY"){
	$PIDsSuper{$ID}{$Data}=1;
	$IsSuper{$Data}=1;
    }elsif(substr($_,0,2) eq "//"){
	$ID="";
    }
}

#Add types to reaction list
foreach my $ID (keys %PIDsTypes){
    foreach my $type (keys %{$PIDsTypes{$ID}}){
	foreach my $rxn (keys %{$ReactionList{$ID}}){
	    $ReactionList{$type}{$rxn}=1;
	}
	traverse_Types_hash($type,\%PIDsTypes,$ID);
    }
}

sub traverse_Types_hash{
    my ($id,$hash,$child)=@_;
    foreach my $key ( grep { exists( $hash->{$_} ) } keys %{$hash->{$id}}){
	foreach my $rxn (keys %{$ReactionList{$child}}){
	    $ReactionList{$key}{$rxn}=1;
	}
	traverse_Types_hash($key,$hash,$child);
    }
}

#Add parents to reaction list
foreach my $ID (keys %PIDsSuper){
    foreach my $super (keys %{$PIDsSuper{$ID}}){
	foreach my $rxn (keys %{$ReactionList{$ID}}){
	    $ReactionList{$super}{$rxn}=1;
	}
	traverse_Supers_hash($super,\%PIDsSuper,$ID);
    }
}

sub traverse_Supers_hash{
    my ($id,$hash,$child)=@_;
    foreach my $key ( grep { exists( $hash->{$_} ) } keys %{$hash->{$id}}){
	foreach my $rxn (keys %{$ReactionList{$child}}){
	    $ReactionList{$key}{$rxn}=1;
	}
	traverse_Supers_hash($key,$hash,$child);
    }
}

my %Rxns_Pwys=();
open(FH, "< ".$Database_Root."/reactions.dat");
($Field,$Data)=("","");
$ID='';
my @Pathways=();
while(<FH>){
    chomp;
    if(substr($_,0,1) eq "#"){next;}
    @temp=split/\s-\s/;
    
    $Field=$temp[0];
    $Data=join("-",@temp[1..$#temp]);
    
    $Field=~s/\s*$//;
    $Data=~s/^\s*//;
    $Data=~s/\s*$//;
    
    if($Field eq "UNIQUE-ID"){
	#07/25
	$ID=$Data;
    }elsif($Field eq "IN-PATHWAY"){
	#some reactions list another reaction as being "in-pathway"
	push(@Pathways,$Data);
    }elsif(substr($_,0,2) eq "//"){
	$Rxns_Pwys{$ID}=() if !exists($Rxns_Pwys{$ID});
	foreach my $p(@Pathways){
	    $Rxns_Pwys{$ID}{$p}=1;
	}
	$ID="";
	undef(@Pathways);
    }
}
close(FH);

#Add pwys for reactions if type or super
foreach my $pwy (sort keys %ReactionList){
    foreach my $rxn (sort keys %{$ReactionList{$pwy}}){
	if(exists($Rxns_Pwys{$rxn})){
	    $Rxns_Pwys{$rxn}{$pwy}=1;
	}
    }
}

#expand pathways used as "reactions"
foreach my $pwy (sort keys %ReactionList){
    foreach my $rxn (sort keys %{$ReactionList{$pwy}}){
	if(!exists($Rxns_Pwys{$rxn})){
	    my $reactions=traverse_PRxn_hash($rxn,\%ReactionList,{});
	    foreach my $rxn2 (keys %$reactions){
		$Rxns_Pwys{$rxn2}{$pwy}=1;
	    }
	}
    }
}

sub traverse_PRxn_hash {
    my ($id, $hash, $reactions) = @_;
    foreach my $key (keys %{$hash->{$id}}){
	$reactions->{$key}=1 if exists($Rxns_Pwys{$key});
	$reactions=traverse_PRxn_hash($key,$hash,$reactions) if !exists($Rxns_Pwys{$key});
    }
    return $reactions;
}

#remove any pathways for any given reaction if it *is* a reaction
foreach my $rxn (sort keys %Rxns_Pwys){
    foreach my $pwy (sort keys %{$Rxns_Pwys{$rxn}}){
	if(exists($Rxns_Pwys{$pwy})){
	    delete($Rxns_Pwys{$rxn}{$pwy});
	}
    }
}

open(OUT, "> MetaCyc_Output/Reactions_Pathways.txt");
foreach my $r (sort keys %Rxns_Pwys){
    foreach my $p(sort keys %{$Rxns_Pwys{$r}}){

	print OUT $r,"\t",$p,"\t";

	print OUT join("|",sort keys %{$PIDsNames{$p}}) if exists($PIDsNames{$p});
	print OUT "\t";

	print OUT join("|",sort keys %{$PIDsTypes{$p}}) if exists($PIDsTypes{$p});
	print OUT "\t";

	if(exists($IsSuper{$p})){
	    print OUT "Super";
	}
	print OUT "\n";
    }
}
close(OUT);

open(OUT, "> ".$Output_Root."MetaCyc_Pathways.tbl");
print OUT "Reaction\tPathway\n";
foreach my $r (sort keys %Rxns_Pwys){
    foreach my $p(sort keys %{$Rxns_Pwys{$r}}){
	next if $p eq "Biosynthesis";
	print OUT $r,"\t",$p,"\t";
	print OUT join("|",sort keys %{$PIDsNames{$p}}) if exists($PIDsNames{$p});
	print OUT "\n";
    }
}
