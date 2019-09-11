#!/usr/bin/env perl
#use PlantSEED::Formulas;
use warnings;
use strict;
my @temp=();

my $Biochemistry = "iJDZ836";
my $Compounds = $Biochemistry."_Compound_Table.txt";
my $Reactions = $Biochemistry."_Reaction_Table.txt";

open(FH, "< $Compounds");
my $header=1;
my %Original_Compounds=();
while(<FH>){
    chomp;
    if($header){$header--;next}
    @temp=split(/\t/,$_);

    my $cpd_cpt = $temp[0];

    #Convert ASCII codes
    $temp[0] =~ s/__43__/+/;
    $temp[0] =~ s/__45__/-/g;
    $temp[0] =~ s/__91__/[/;
    $temp[0] =~ s/__93__/]/;

    #Clean up identifier
    $temp[0] =~ s/^M_+//;

    #Remove Compartment
    my $cpd = $temp[0];

    my $cpt="c";
    if($cpd =~ /\[/){
	$cpd =~ s/([\w+-]+)\[([\w-]+?)\]$/$1/;
	if(defined($2)){
	    $cpt = $2;
	}
    }elsif($cpd =~ /_bm$/){
	$cpd =~ s/([\w+-]+)_(\w+?)$/$1/;
	if(defined($2)){
	    $cpt = $2;
	}
    }
    $Original_Compounds{$cpd_cpt}={'ID'=>$cpd,
				   'NAMES'=>$temp[1],
				   'COMPARTMENT'=>$cpt};
}
close(FH);

open(FH, "< $Reactions");
$header=1;
my %Original_Reactions=();
my %Cpds_in_Rxns=();
while(<FH>){
    chomp;
    if($header){$header--;next}

    #Convert ASCII codes
    $_ =~ s/__43__/+/g;
    $_ =~ s/__45__/-/g;
    $_ =~ s/__91__/[/g;
    $_ =~ s/__93__/]/g;

    @temp=split(/\t/,$_);

    my $rxn = $temp[0];

    #Clean up identifier
    $rxn =~ s/^R//;

    #Go through reactants
    my ($rev,$reactants,$products)=@temp[2..5];

    #Skipping all boundary exchange reactions
    next if !$products;
    
    my @reactants = split(/;/,$reactants);

    my @eqn=();
    my %cpts = (); #got to double-check
    my $cpt_count = 0;
    foreach my $rct (@reactants){
	$rct =~ s/^M_//;
	my ($cpd,$cpt,$coeff)=($rct,"c","0");
	$rct =~ s/\[([\d\.]+)\]$//;
	$coeff=$1;
	$cpd=$rct;
	if($rct =~ /\[/){
	    $cpd =~ s/([\w+-]+)\[([\w-]+?)\]$/$1/;
	    if(defined($2)){
		$cpt = $2;
	    }
	}elsif($rct =~ /_bm$/){
	    $cpd =~ s/([\w+-]+)_(\w+?)$/$1/;
	    if(defined($2)){
		$cpt = $2;
	    }
	}

	$Cpds_in_Rxns{$cpd}=1;

	if(defined($coeff) && $coeff != 1){
	    $coeff="(".$coeff.")";
	}else{
	    $coeff=undef($coeff);
	}

	if(!exists($cpts{$cpt})){
	    $cpts{$cpt}=$cpt_count;
	    $cpt_count++;
	}

	$cpt = "[".$cpts{$cpt}."]";
	
	my $rct_str = $cpd.$cpt;
	if(defined($coeff)){
	    $rct_str = $coeff." ".$rct_str;
	}

	push(@eqn,$rct_str);
    }

    my $reversibility = "<=>";
    if($rev ne "False"){
	$reversibility = "=>";
    }
    push(@eqn,$reversibility);

    my @products = split(/;/,$products);
    foreach my $pdt (@products){

	$pdt =~ s/^M_//;
	my ($cpd,$cpt,$coeff)=($pdt,"c","0");
	$pdt =~ s/\[([\d\.]+)\]$//;
	$coeff=$1;
	$cpd=$pdt;
	if($pdt =~ /\[/){
	    $cpd =~ s/([\w+-]+)\[([\w-]+?)\]$/$1/;
	    if(defined($2)){
		$cpt = $2;
	    }
	}elsif($pdt =~ /_bm$/){
	    $cpd =~ s/([\w+-]+)_(\w+?)$/$1/;
	    if(defined($2)){
		$cpt = $2;
	    }
	}

	$Cpds_in_Rxns{$cpd}=1;

	if(defined($coeff) && $coeff != 1){
	    $coeff="(".$coeff.")";
	}else{
	    $coeff=undef($coeff);
	}

	if(!exists($cpts{$cpt})){
	    $cpts{$cpt}=$cpt_count;
	    $cpt_count++;
	}

	$cpt = "[".$cpts{$cpt}."]";

	my $pdt_str = $cpd.$cpt;
	if(defined($coeff)){
	    $pdt_str = $coeff." ".$pdt_str;
	}
	push(@eqn,$pdt_str);
    }
    my $eqn_str = join(" ",@eqn);

    $Original_Reactions{$rxn}={'ID'=>$rxn,
			       'NAMES'=>$temp[1],
			       'EQUATION'=>$eqn_str};
}
close(FH);

my $filestub = $Compounds;
$filestub =~ s/_Compound_Table\.txt$//;

open(OUT, "> ".$filestub."_Compounds.tbl");
my @Headers=("ID","NAMES","COMPARTMENT");
print OUT join("\t",@Headers),"\n";
foreach my $id (sort keys %Original_Compounds){
    next if !exists($Cpds_in_Rxns{$Original_Compounds{$id}{'ID'}});

    foreach my $h (@Headers){
	print OUT $Original_Compounds{$id}{$h};
	print OUT "\t" unless $h eq $Headers[$#Headers];
    }
    print OUT "\n";
}
close(OUT);

open(OUT, "> ".$filestub."_Reactions.tbl");
@Headers=("ID","NAMES","EQUATION");
print OUT join("\t",@Headers),"\n";
foreach my $id (sort keys %Original_Reactions){
    foreach my $h (@Headers){
	print OUT $Original_Reactions{$id}{$h};
	print OUT "\t" unless $h eq $Headers[$#Headers];
    }
    print OUT "\n";
}
close(OUT);
