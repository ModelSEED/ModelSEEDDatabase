package Reactions;
use strict;
use warnings;

sub balance_reaction {
    my $cpd_hash = shift;
    if(scalar(keys %$cpd_hash)==0){
	return "EMPTY";
    }

    my @status = ("OK");

    #Code aimed at finding duplicates
    my %Rgts_Hash=();
    foreach my $cpd (keys %$cpd_hash){
	$Rgts_Hash{$cpd."_".$cpd_hash->{$cpd}{compartment}}++;
    }
    if(scalar( grep { $Rgts_Hash{$_} > 1 } keys %Rgts_Hash)>0){
	#Here reagents removed
    }

    foreach my $cpd (keys %$cpd_hash){
	$Rgts_Hash{$cpd."_".$cpd_hash->{$cpd}{compartment}}++;
    }

    my $protonCompHash=();
    my $compHash=();

    #check for reactions with duplicate reagents in different compartments
    #these are not removed, but they balance out regardless of formula
    my %Cpds_Hash = ();
    foreach my $rgt (sort keys %Rgts_Hash){
	my ($cpd,$cpt) = split(/_/,$rgt);

	#Check for protons/water
	$protonCompHash->{$cpt}=1 if checkForProton($cpd);

	#Store compartments
	$compHash->{$cpt}=1;

	#balance out reagents regardless of compartment
	$Cpds_Hash{$cpd}+=$cpd_hash->{$cpd}{coefficient};
    }

    my $atomHash={};
    my $netCharge=0;
    foreach my $cpd_id ( grep { $Cpds_Hash{$_} !=0 } keys %Cpds_Hash){

	#Problems are: compounds with noformula/null and polymers
	#Latest KEGG formulas for polymers contain brackets and 'n', older ones contain '*'
	my $cpdatoms = Formulas::parse($cpd_hash->{$cpd_id}{'formula'});
	return "CPDFORMERROR" if exists($cpdatoms->{error});
	
	#Add up charge
	$netCharge += $cpd_hash->{$cpd_id}{'charge'}*$cpd_hash->{$cpd_id}{'coefficient'};
	
	#Add up atoms
	foreach my $atom (keys(%{$cpdatoms})) {
	    if (!defined($atomHash->{$atom})) {
		$atomHash->{$atom} = 0;
	    }
	    $atomHash->{$atom} += $Cpds_Hash{$cpd_id}*$cpdatoms->{$atom};
	}
    }

    #Add protons if missing
    if (!defined($atomHash->{H})) {
	$atomHash->{H} = 0;
    }

    my $imbalancedAtoms = {};
    foreach my $atom (keys(%{$atomHash})) { 
	if ($atomHash->{$atom} > 0.00000001 || $atomHash->{$atom} < -0.00000001) {
	    $imbalancedAtoms->{$atom}=$atomHash->{$atom};
	}
    }

    if (join("",keys %$imbalancedAtoms) eq "H") {
	if(scalar(keys %$protonCompHash)==0){
	    #defaults to 'highest' compartment
	    my $first_comp = (sort keys %$compHash)[0];

	    #Add proton
	    $cpd_hash->{'cpd00067'}={'formula'=> 'H',
				     'charge' => 1,
				     'coefficient' => -1*$imbalancedAtoms->{"H"},
				     'compartment' => $first_comp};

	}elsif(scalar(keys %$protonCompHash)>0){
	    #defaults to 'highest' compartment
	    my $first_comp = (sort keys %$protonCompHash)[0];

	    my $coeff = $cpd_hash->{'cpd00067'}{'coefficient'};
	    $coeff = $coeff + (-1*$imbalancedAtoms->{"H"});
	    $cpd_hash->{'cpd00067'}{'coefficient'}=$coeff;
	}
	
	#adjust charge based on proton balancing
	$netCharge += -1*$imbalancedAtoms->{"H"};

	#balance out protons
	$atomHash->{H} = 0;
	delete($imbalancedAtoms->{H});
	push(@status,"HB");
    }

    #fix tiny coefficients
    foreach my $atom (keys %{$imbalancedAtoms}){
	if($atomHash->{$atom} > -0.00000001 && $atomHash->{$atom} < 0.00000001) {
	    $atomHash->{$atom}=0;
	    delete($imbalancedAtoms->{$atom});
	}
    }

    #report imbalance
    if(scalar(keys %{$imbalancedAtoms})>0){
	$status[0] = "MI:".join("/", map { $_.":".$atomHash->{$_} } sort keys %{$imbalancedAtoms});	
    }

    #fix tiny charge
    if($netCharge > -0.00000001 && $netCharge < 0.00000001) {
	$netCharge=0;
    }
    
    if ($netCharge != 0) {
	push(@status,"CI:".$netCharge);
	shift(@status) if $status[0] eq "OK";
    }

    return join("|",@status);
}

sub checkForProton {
    my $id = shift;
    return $id =~ /^(cpd00067|C00080|PROTON)$/;
}
1;
