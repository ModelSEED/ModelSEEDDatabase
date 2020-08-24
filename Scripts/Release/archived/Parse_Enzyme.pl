#!/usr/bin/env perl
use warnings;
use strict;

#ID   1.1.1.1
#DE   Alcohol dehydrogenase.
#AN   Aldehyde reductase.
#CA   (1) A primary alcohol + NAD(+) = an aldehyde + NADH.
#CA   (2) A secondary alcohol + NAD(+) = a ketone + NADH.

my %Enzymes=();
my ($ID,$Names,$Equations)=(undef,[],[]);

open(FH, "< enzyme.dat");
while(<FH>){
    chomp;

    $_ =~ s/(^[\/\w]{2})//;
    my $line_type = $1;
    my $line = $_;

    #Skip comments
    next if $line_type eq "CC";

    #strip whitespace, for all keys lines, it is always three spaces.
    $line =~ s/^\s+//;

    #there is never trailing whitespace, so this is redundant
    $line =~ s/\s+$//;

    if($line_type eq "ID"){
	$ID=$line;
    }
    if($line_type eq "DE" or $line_type eq "AN"){
	if(scalar(@$Names)>0 && $Names->[-1] !~ /\.$/){
	    $Names->[-1].=$line;
	}else{
	    push(@$Names,$line);
	}
    }
    if($line_type eq "CA"){
	if(scalar(@$Equations)>0 && $Equations->[-1] !~ /\.$/){
	    if($Equations->[-1] =~ /[\+=]$/ || $line =~ /^[\+=]/){
		$Equations->[-1].=" ".$line;
	    }else{
		$Equations->[-1].=$line;
	    }
	}else{
	    push(@$Equations,$line);
	}
    }
    if($line_type eq "//"){
	if($ID && scalar(@$Equations)>0){
	    $Enzymes{$ID}{"names"}=[@$Names];
	    $Enzymes{$ID}{"equations"}=[@$Equations];
	}
	($ID,$Names,$Equations)=(undef,[],[]);
    }
}
close(FH);

open(OUT, "> Parsed_Enzyme_Equations.txt");
foreach my $id (sort keys %Enzymes){
    foreach my $eqn (@{$Enzymes{$id}{'equations'}}){
	print OUT $id,"\t",$eqn,"\n";
    }
}
