#!/usr/bin/env perl
use warnings;
use strict;
use Formulas;
my @temp=();

my $File = $ARGV[0];
die("No file or $File isn't found") if !$File || !-f $File;

#recover modifications that have already been applied to avoid duplicating records
open(FH, "< ../Biochemistry/compounds.master.mods");
my %Touched=();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $Touched{$temp[0]}=1;
}
close(FH);

open(FH, "< $File");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next if exists($Touched{$temp[0]});

    my $formula1 = Formulas::merge_formula($temp[1]);
    my $formula2 = Formulas::merge_formula($temp[2]);

    #empty formulas
    if($temp[1] eq "" && ( $temp[2] eq "unknown" || $temp[2] eq "noformula") ){
#	print $temp[0],"\tcurated\tformula\tnull\n";
	print "EMPTY\t".$_."\n";
	next;
    }

    if($temp[1] eq "" || $temp[2] eq "unknown" || $temp[2] eq "noformula"){
	#formula missing altogether from one db
	if($temp[1] eq ""){
#	    print $temp[0],"\tplantdefault\tformula\t",$formula2,"\n";
	}elsif($temp[2] eq "unknown" || $temp[2] eq "noformula"){
#	    print $temp[0],"\tdefault\tformula\t",$formula1,"\n";
	}
	
	print "MISSING\t".$_."\n";
	next;
    }

    if($formula1 eq $formula2){
#	print "MATCH\t".$_."\n";
	next;
    }

    my $formula1hash = Formulas::parse($formula1);
    my $formula2hash = Formulas::parse($formula2);

    if(scalar( grep { !exists($formula2hash->{$_}) } keys %$formula1hash)>0 || scalar( grep { !exists($formula1hash->{$_}) } keys %$formula2hash)>0){
	if( ( scalar( grep { !exists($formula2hash->{$_}) } keys %$formula1hash )>0 && ( grep { !exists($formula2hash->{$_}) } keys %$formula1hash )[0] eq "H" ) || 
	    ( scalar( grep { !exists($formula1hash->{$_}) } keys %$formula2hash )>0 && ( grep { !exists($formula1hash->{$_}) } keys %$formula2hash )[0] eq "H" ) ){
#	    print $temp[0],"\tplantdefault\tformula\t",$formula2,"\n";
	}

	#missing elements
	print "ELEMENTS\t".$temp[0],"\t",join("|", sort grep { !exists($formula2hash->{$_}) } keys %$formula1hash)."\t".join("|", sort grep { !exists($formula1hash->{$_}) } keys %$formula2hash)."\n";
	next;
    }

    my %diffhash = map { $_ => $formula1hash->{$_}-$formula2hash->{$_} } keys %$formula1hash;
    if(scalar( grep { $diffhash{$_} != 0 } keys %diffhash)>0){
	if(scalar(grep {$diffhash{$_}!=0} keys %diffhash)==1 && (grep {$diffhash{$_}!=0} keys %diffhash)[0] eq "H"){
#	    print $temp[0],"\tplantdefault\tformula\t",$formula2,"\n";
	}


	print "NUMBERS\t".$temp[0],"\t",join("|", map { $_.":".$diffhash{$_} } grep { $diffhash{$_} !=0 } keys %diffhash),"\n";
    }
}
