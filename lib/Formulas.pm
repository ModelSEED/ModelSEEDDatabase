package Formulas;
use strict;
use warnings;

sub parse{
    my $string=shift;
    return {} if !$string;
    return {error => "No formula"} if $string =~ /^(noformula|null)$/;
    my @atoms = ($string =~ /(\D[a-z]?\d*)/g);
    my %atomKey=();
    foreach my $a (@atoms){
        my @stoich=($a =~ /(\D[a-z]?)(\d*)/);
        $stoich[1]=1 if $stoich[1] eq '';
        $atomKey{$stoich[0]}=$stoich[1];
    }
    return \%atomKey;
}

sub merge_formula{
    my $formula=shift;
    my $remove=shift;

    #polymers
    return $formula if $formula =~ /\*\d?$/;

    #remove spaces
    $formula =~ s/\s//g;

    while($formula =~ /\((.*?)\)([\w\*]*)/g){
	my $bracket_multi=$2;
	$bracket_multi=1 if !$bracket_multi;
	if($bracket_multi !~ /^\d+$/){
	    return $formula;
	}
    }

    my %global_atoms=();    
    while($formula =~ /([\w\s\.]*)[\(\)]?(\d?)/g){
	if($1){
	    my $string=$1;
	    my $bracket_multi=1;
	    $bracket_multi=$2 if $2;

	    foreach my $entity ( grep { $_ ne "" } split(/\./,$string)){

		#get next multiple
		my $string_multi=1;
		if($entity =~ s/^(\d+)(.*)$/$2/){
		    $string_multi=$1;
		}

		my $atoms=parse($entity);
		foreach my $a(keys %$atoms){
		    $global_atoms{$a}+=($atoms->{$a}*$bracket_multi*$string_multi);
		}
	    }
	}
    }
    if($remove){
	foreach my $element (@$remove){
	    delete($global_atoms{$element});
	}
    }
    return hill_sort_formula(\%global_atoms);
}

sub hill_sort_formula{
    my $atom_hash=shift;

    my @atom_order=();
    my $exclude='\*';

    #Hill sorted formula
    #http://en.wikipedia.org/wiki/Hill_system

    if(exists($atom_hash->{'C'})){
	push(@atom_order,'C');
	push(@atom_order,'H') if exists($atom_hash->{'H'});
	$exclude.='CH';
    }
    push(@atom_order, grep { $_ !~ /^[$exclude]$/ } sort keys %$atom_hash );
    push(@atom_order,'*') if exists($atom_hash->{'*'});

    my $new_formula="";
    foreach my $a (@atom_order){
	$atom_hash->{$a}="" if $atom_hash->{$a} && $atom_hash->{$a} == 1;
	$new_formula.=$a.$atom_hash->{$a};
    }

    return $new_formula;
}
1;
