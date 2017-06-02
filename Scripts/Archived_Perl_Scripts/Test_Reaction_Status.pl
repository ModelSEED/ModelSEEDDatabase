#!/usr/bin/env perl
use warnings;
use strict;
use Formulas;
my @temp=();
my $header = 1;

my $Test_Rxn = $ARGV[0];
exit if !$Test_Rxn;

my $db="master";
$db = $ARGV[1] if $ARGV[1];

#######################################################
#Create Empty Biochemistry
#######################################################
use Bio::KBase::ObjectAPI::KBaseBiochem::Biochemistry;
my $Bio_Obj=Bio::KBase::ObjectAPI::KBaseBiochem::Biochemistry->new({id=>"Temp"});
$Bio_Obj->_reference("~");

#######################################################
#Load Master Compounds
#######################################################
open(FH, "< ../Biochemistry/compounds.".$db.".tsv");
my @headers = split(/\t/,<FH>);
chomp($headers[$#headers]);
my %Required_Headers=(id=>'id',charge=>'defaultCharge',formula=>'formula',name=>'name');
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);

    my %Cpd_Hash = map { $Required_Headers{$headers[$_]} => $temp[$_] } grep { exists($Required_Headers{$headers[$_]}) } (0..$#temp);

    #Reverting to original default values
    $Cpd_Hash{defaultCharge}=10000000 if $Cpd_Hash{defaultCharge} eq "null";
    $Cpd_Hash{formula}="noformula" if $Cpd_Hash{formula} eq "null";

    my $Obj = $Bio_Obj->add("compounds",\%Cpd_Hash);
}
close(FH);

#######################################################
#Iterate Master Reactions
#######################################################
open(FH, "< ../Biochemistry/reactions.".$db.".tsv");
@headers = split(/\t/,<FH>);
chomp($headers[$#headers]);
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next unless $temp[0] eq $Test_Rxn;

    my %Rxn_Hash = map { $headers[$_] => $temp[$_] } (0..$#temp);
    my $Rxn_Obj = $Bio_Obj->add("reactions",{id=>$temp[0]});

    foreach my $cpd_array (split(/\;/,$Rxn_Hash{stoichiometry})){
	my ($coef,$cpd,$cmpt) = split(/:/,$cpd_array);

	$cmpt = "0";

	#Check for compartment
	my $cmptObj = $Bio_Obj->getObject("compartments", $cmpt);
	unless(defined($cmptObj)) {
	    $cmptObj = $Bio_Obj->add("compartments",{id => $cmpt,hierarchy=>3});
	}

	#Add reagent
	$Rxn_Obj->add("reagents",{compound_ref => "~/compounds/id/".$cpd,
				  compartment_ref => "~/compartments/id/".$cmpt,
				  coefficient => $coef,
				  isCofactor => 0});

    }

    print $Rxn_Hash{equation}."\t".$Rxn_Hash{stoichiometry}."\t".$Rxn_Hash{status}."\n";

    $Rxn_Obj->status($Rxn_Hash{status});
    $Rxn_Obj->checkReactionMassChargeBalance({rebalanceProtons=>1,rebalanceWater=>0,saveStatus=>1});

    print $Rxn_Obj->equation()."\t".$Rxn_Obj->stoichiometry()."\t".$Rxn_Obj->status()."\n";
}
close(FH);
