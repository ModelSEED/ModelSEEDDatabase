#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header=1;

open(FH, "< Source_Classifiers.txt");
$header=1;
my %Source_Classes=();
while(<FH>){
    chomp;
    if($header){$header--;next}
    my ($source_id,$source_class)=split(/\t/,$_);
    $Source_Classes{$source_id}=$source_class;
}
close(FH);

#Retrieve necessary aliases
my %Sources=();
foreach my $type ('Compounds','Reactions'){
    foreach my $stub ('_Aliases.tsv'){
	my $file = $type.$stub;
	open(FH, "< $file");
	$header=1;
	while(<FH>){
	    chomp;
	    if($header){$header--;next}
	    @temp=split(/\t/,$_);

	    if($temp[0] =~ /(cpd|rxn)/){
		foreach my $id (split(/\|/,$temp[0])){
		    $Sources{$id}{$temp[3]}=1;
		}
	    }elsif($temp[1] =~ /(cpd|rxn)/){
		foreach my $id (split(/\|/,$temp[1])){
		    $Sources{$id}{$temp[3]}=1;
		}
	    }else{
		print $_,"\n";
	    }
	}
	close(FH)
    }
}

my $Type="compounds";
open(FH, "< $ENV{SEAVER_PROJECT}ModelSEEDDatabase/Biochemistry/".$Type.".tsv");
my %Biochemistry=();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $Biochemistry{$temp[0]}=$_;
}
close(FH);

#Classify
foreach my $id (sort keys %Biochemistry){
    next if $id eq "id";
    my ($primary,$secondary,$published)=(undef,undef,undef);
    if(exists($Sources{$id})){
	foreach my $source (sort keys %{$Sources{$id}}){
	    if($Source_Classes{$source} eq "Primary Database"){
		$primary=$Source_Classes{$source};
	    }
	    if($Source_Classes{$source} eq "Secondary Database"){
		$secondary=$Source_Classes{$source};
	    }
	    if($Source_Classes{$source} eq "Published Model"){
		$published=$Source_Classes{$source};
	    }
	}
    }

    my $source=undef;
    if(defined($primary)){
	$source=$primary;
    }elsif(defined($secondary)){
	$source=$secondary;
    }elsif(defined($published)){
	$source=$published;
    }else{
	$source="Orphan";
    }
    $Biochemistry{$id}.="\t".$source;
}

$Biochemistry{id}.="\tsource";

open(OUT, "> $ENV{SEAVER_PROJECT}ModelSEEDDatabase/Biochemistry/".$Type.".old");
print OUT $Biochemistry{id}."\n";
foreach my $id (sort keys %Biochemistry){
    next if $id eq "id";
    print OUT $Biochemistry{$id},"\n";
}
close(OUT);
