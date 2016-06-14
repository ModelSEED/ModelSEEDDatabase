package InChIs;
use strict;
use warnings;
use Formulas;

my $ProjectPath=$ENV{SEAVER_PROJECT}."Structures/";
my @InChI_Layers=('c','h','p','q','b','t','m','s');
my @temp=();

sub layers{
    return \@InChI_Layers;
}

sub get_All_InChIs{
    my $args=shift;

    $args->{type}="Original";
    my $InChIs=get_InChIs($args);

    $args->{type}="Charged";
    my $inchis = get_InChIs($args);

    foreach my $ma (keys %$inchis){
	$InChIs->{$ma}=$inchis->{$ma};
    }

    return $InChIs;
}

sub get_InChIs{
    my $self = shift if $_[0] eq "InChIs";
    my $args=shift;
    my %InChIs=();

    my $file=$ProjectPath."InChI/".$args->{source}."/".$args->{type}."Strings.txt";
    if(!-f $file){
	print "$file doesn't exist\n";
	return;
    }

    open(FH, "< $file");
    while(<FH>){
	chomp;
	@temp=split(/\t/);
	$InChIs{$temp[0]}=$temp[1];
    }

    return \%InChIs;
}

sub parse{
    my $inchi=shift;
    my $args=shift;

    my @layers=split /\//, $inchi;
    shift @layers; #remove InChI
    my $formula = shift @layers;
 
    my @components=();
    foreach my $component(split(/\./,$formula)){
	$component=Formulas::merge_formula($component);
	push(@components,$component);
    }
    $formula=join(".",@components);

    my %layers=();
    foreach my $l (@InChI_Layers){
	$layers{$l}='';
    }

    foreach my $l (@layers){
        $l =~ /(.)(.*)/;
        my ($layer_id, $layer_data) = ($1,$2);
        $layers{$layer_id}=$layer_data;
    }

    #special fix for C00080 which doesn't have a formula
    #but a p+1 layer (yea, its just a proton)
    if($inchi =~ /^InChI=1S\/p([-+]\d*)/){
        $layers{'p'}=$1;
        $formula="";
    }

    $formula=Formulas::merge_formula($formula) if exists($args->{merge_formula});

    return ($formula,\%layers);
}

sub build{
    my ($formula,$layers,$args)=@_;

    $formula=Formulas::merge_formula($formula) if exists($args->{merge_formula});
    
    #re-build inchi string
    #the "search" version has no "p" and "q" layers

    my $string="InChI=1S";
    
    #special case of single proton not having a formula
    $string.="/".$formula if $formula;

    foreach my $l( grep { $layers->{$_} ne '' } @{InChIs::layers()}){
	$string.="/".$l.$layers->{$l} unless exists($args->{remove}{$l});
    }

    return "" if $string eq "InChI=1S"; #proton case again but possibly other erroneous ones
    return $string;
}

sub charge{
    my ($q_string,$p_string)=@_;

    my $charge=0;
    foreach my $q ( grep { $_ ne "" } split(/;/, $q_string)){
	my $multiple=1;
	if($q =~ s/^(\d+)\*(.+)$/$2/){
	    $multiple=$1;
	}
	$charge+=($q*$multiple);
    }

    foreach my $p ( grep { $_ ne "" } split(/;/, $p_string)){
	my $multiple=1;
	#NB: at time of press (12/2/2011)
	#No /p sublayer does multiple components
	#therefore multiple code doesn't apply
	if($p =~ s/^(\d+)\*(.+)$/$2/){
	    $multiple=$1;
	}
	$charge+=($p*$multiple);
    }
    return $charge;
}

sub adjust_protons{
    my ($formula,$protons)=@_;

    #special case of protons
    if($formula eq ""){
	return "H";
    }

    return $formula if !$protons;

    my @components=split(/\./,$formula);

    #NB: at time of press (3/19/2015)
    #All protons are accounted for in a single component
    #Mostly the first component, but rarely the second component where the first does not contain protons
    #Proton adjustment within a single component never goes below zero

    for(my $i=0;$i<scalar(@components);$i++){
	my $atoms=Formulas::parse($components[$i]);
	if($atoms->{"H"}){
	    $atoms->{"H"}+=$protons;

	    if($atoms->{"H"} < 0){
		print STDERR "ERROR: Too Many Protons adjusted in formula!\n";
	    }

	    delete($atoms->{"H"}) if $atoms->{"H"}==0;
	    $components[$i]=Formulas::hill_sort_formula($atoms);
	    last;
	}
    }

    return join(".",@components);
}
1;

__END__

#$AuxInfoString=$_ if $_ =~ /^AuxInfo=/;
#if($AuxInfoString =~ /CRV/){
#    $Radicals{$KEGG}=1;
#}
