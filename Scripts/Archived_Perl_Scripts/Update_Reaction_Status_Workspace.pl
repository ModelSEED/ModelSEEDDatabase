#!/usr/bin/env perl
use warnings;
use strict;
use Formulas;
use Getopt::Long::Descriptive;
my @temp=();
my $header = 1;

my ($opt, $usage) = describe_options("%c %o ",
	[ "conserve_hb", "conserver hb in reaction status" ],
	[ "help|h", "print usage message and exit" ]
);
print($usage->text), exit if $opt->help;

#######################################################
#Create Empty Biochemistry
#######################################################
use Bio::KBase::ObjectAPI::KBaseBiochem::Biochemistry;
my $Bio_Obj=Bio::KBase::ObjectAPI::KBaseBiochem::Biochemistry->new({id=>"Temp"});
$Bio_Obj->_reference("~");

#######################################################
#Load Master Compounds
#######################################################
open(FH, "< ../Biochemistry/compounds.master.tsv");
my @headers = split(/\t/,<FH>);
chomp($headers[$#headers]);
my %Required_Headers=(id=>'id',charge=>'defaultCharge',formula=>'formula',name=>'name');
my %Cpds=();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);

    my %Cpd_Hash = map { $Required_Headers{$headers[$_]} => $temp[$_] } grep { exists($Required_Headers{$headers[$_]}) } (0..$#temp);

    #Reverting to original default values
    $Cpd_Hash{defaultCharge}=10000000 if $Cpd_Hash{defaultCharge} eq "null";
    $Cpd_Hash{formula}="noformula" if $Cpd_Hash{formula} eq "null";

    my $Obj = $Bio_Obj->add("compounds",\%Cpd_Hash);
    $Cpds{$temp[0]}=$Obj;
}
close(FH);

#Retrieve prioritized reactions
open(FH, "< ../Biochemistry/Workspaces/KBaseTemplateModels.rxn");
my %PriRxns=();
while(<FH>){
    chomp;
    @temp=split(/\t/);
    $PriRxns{$temp[0]}=1;
}
close(FH);

#######################################################
#Iterate Master Reactions
#######################################################
open(FH, "< ../Biochemistry/reactions.master.tsv");
@headers = split(/\t/,<FH>);
chomp($headers[$#headers]);
my %Rxns=();
open(OUT, "> Status_Diffs.txt");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
#    next unless $temp[0] eq "rxn05701";

    my %Rxn_Hash = map { $headers[$_] => $temp[$_] } (0..$#temp);
    my $Rxn_Obj = $Bio_Obj->add("reactions",{id=>$temp[0]});

    foreach my $cpd_array (split(/\;/,$Rxn_Hash{stoichiometry})){
	my ($coef,$cpd,$cmpt) = split(/:/,$cpd_array);

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

    if ($opt->conserve_hb){
	$Rxn_Obj->status($Rxn_Hash{status});
    }
    $Rxn_Obj->checkReactionMassChargeBalance({rebalanceProtons=>1,rebalanceWater=>0,saveStatus=>1});

    print OUT exists($PriRxns{$Rxn_Hash{id}})."\t".$Rxn_Hash{id}."\t".$Rxn_Hash{status}."\t".$Rxn_Obj->status()."\n";

    $Rxn_Hash{equation}=$Rxn_Obj->genEquation();
    $Rxn_Hash{code}=$Rxn_Obj->genCode();
    $Rxn_Hash{definition}=$Rxn_Obj->genDefinition();
    $Rxn_Hash{stoichiometry}=$Rxn_Obj->genStoichiometry();
    $Rxn_Hash{status}=$Rxn_Obj->status();

    $Rxns{$temp[0]}=join("\t", map { $Rxn_Hash{$_} } @headers);
}
close(FH);
close(OUT);

open(OUT, "> ../Biochemistry/reactions.master.tsv");
#open(OUT, "> tmp");
print OUT join("\t",@headers)."\n";
foreach my $rxn ( sort { $a cmp $b } keys %Rxns){
    print OUT $Rxns{$rxn},"\n";
}
close(OUT);

__END__

#Old Status Checking Code


#Retrieve prioritized reactions
open(FH, "< ../Biochemistry/Workspaces/KBaseTemplateModels.rxn");
my %PriRxns=();
while(<FH>){
    chomp;
    @temp=split(/\t/);
    $PriRxns{$temp[0]}=1;
}
close(FH);

#Retrieve list of modified compounds
open(FH, "< ../Biochemistry/compounds.master.mods");
my %Compounds=();
while(<FH>){
    chomp;
    @temp=split(/\t/);
    $Compounds{$temp[0]}{$temp[2]}=$temp[3] if $temp[2] eq "charge" || $temp[2] eq "formula";
}
close(FH);

use Bio::KBase::fbaModelServices::ScriptHelpers qw( getToken );
my $AToken = getToken();

use Bio::KBase::fbaModelServices::Impl;
my $FBAImpl = Bio::KBase::fbaModelServices::Impl->new({'fbajobcache' => "/homes/seaver/Projects/KBase_Scripts/FBA_Scripts/JobsCache",
						       'jobserver-url' => "http://kbase.us/services/workspace",
						       'fbajobdir' => "/tmp/fbajobs",
						       'mfatoolkitbin' => "/vol/model-prod/kbase/MFAToolkit/bin/mfatoolkit",
#						       'mfatoolkitbin' => "/homes/seaver/Software/MFAToolkit/bin/mfatoolkit",
						       'probanno-url' => "http://140.221.85.86:7073/",
						       'mssserver-url' => "http://bio-data-1.mcs.anl.gov/services/ms_fba",
						       'accounttype' => "kbase",
						       'workspace-url' => "http://kbase.us/services/ws",
						       'defaultJobState' => "queued",
						       'gaserver-url' => "http://kbase.us/services/genome_annotation",
						       'idserver-url' => "http://kbase.us/services/idserver"});
$FBAImpl->_setContext(undef,{auth=>$AToken});
my $Store = $FBAImpl->_KBaseStore();

my $dObj = $FBAImpl->_get_msobject("Biochemistry", 'kbase', 'default');
my $pdObj = $FBAImpl->_get_msobject("Biochemistry", 'kbase', 'plantdefault');

#foreach my $rxn ( grep { scalar( @{$_->getAliases("ModelSEED")} ) > 1 } @{$pdObj->reactions()}){
#    print $rxn,"\t",join("|", sort @{$rxn->getAliases("ModelSEED")}),"\n";
#}

#Differing reaction codes are found where compound meaning (and their id) had been changed
#For this present time, we're defaulting to the reaction in default, and making plantdefault cross-compatible
#So I'm collecting the ones that differ, for curation in the near future, but ignoring them for now
my %Differing_Reactions=();
foreach my $rxn (sort keys %PriRxns){
    my $dRxn = $dObj->getObject("reactions",$rxn);
    my $pdRxn = $pdObj->getObject("reactions",$rxn);

    if($dRxn && $pdRxn){
	if($dRxn->equationCompFreeCode() ne $pdRxn->equationCompFreeCode() && $dRxn->equationCompFreeCode() ne $pdRxn->revEquationCompFreeCode()){
	    $Differing_Reactions{$rxn}=1;
	}
    }
}

my $cpd_loc_rxn_hash = $dObj->compound_reaction_hash();

my %Reactions=();
foreach my $cpd ( grep { exists($Compounds{$_}) } keys %$cpd_loc_rxn_hash ){
    foreach my $loc (keys %{$cpd_loc_rxn_hash->{$cpd}}){
	foreach my $rxn ( grep { exists ($PriRxns{$_}) } keys %{$cpd_loc_rxn_hash->{$cpd}{$loc}}){
	    next if exists($Differing_Reactions{$rxn});
	    $Reactions{$rxn}{$cpd}=1;
	}
    }
}

foreach my $rxn (keys %Reactions){
    my $dRxn = $dObj->getObject("reactions",$rxn);
    my $pdRxn = $pdObj->getObject("reactions",$rxn);

    #In very few cases, due to changes in compounds, reactions were merged, and the reaction id in default disappeared from plantdefault
    next if !$pdRxn;

    #Here Im double-checking whether the formula or charge need to be updated before balancing
    foreach my $rgt ( grep { exists($Compounds{$_->compound()->id()}) } @{$dRxn->reagents()} ){
	my $Cpd = $Compounds{$rgt->compound()->id()};

	if(exists($Cpd->{formula})){
	    if($Cpd->{formula} ne $rgt->compound()->formula()){
		$rgt->compound()->formula($Cpd->{formula});
	    }
	}

	if(exists($Cpd->{charge})){
	    if($Cpd->{charge} ne $rgt->compound()->defaultCharge()){
		$rgt->compound()->defaultCharge($Cpd->{charge});
	    }
	}
    }

    #And again for the reaction in plantdefault, remember they have the same code
    foreach my $rgt ( grep { exists($Compounds{$_->compound()->id()}) } @{$pdRxn->reagents()} ){
	my $Cpd = $Compounds{$rgt->compound()->id()};

	if(exists($Cpd->{formula})){
	    if($Cpd->{formula} ne $rgt->compound()->formula()){
		$rgt->compound()->formula($Cpd->{formula});
	    }
	}

	if(exists($Cpd->{charge})){
	    if($Cpd->{charge} ne $rgt->compound()->defaultCharge()){
		$rgt->compound()->defaultCharge($Cpd->{charge});
	    }
	}
    }

    #Rebalancing both reactions, and testing the change in status
    my $olddStatus= $dRxn->status();
    $dRxn->checkReactionMassChargeBalance({rebalanceProtons=>1,rebalanceWater=>0,saveStatus=>1});
    my $newdStatus = $dRxn->status();

    my $oldpdStatus= $pdRxn->status();
    $pdRxn->checkReactionMassChargeBalance({rebalanceProtons=>1,rebalanceWater=>0,saveStatus=>1});
    my $newpdStatus = $pdRxn->status();

    #Special case when all four is OK, and we default to the new status of the plantdefault reaction as it's the most recent
    if( scalar( grep { $_ =~ /^OK/ } ( $olddStatus,$newdStatus,$oldpdStatus, $newpdStatus )) == 4 ){
	
	#Currently, if proton balancing was done that changes the charge-balance of the reaction
	#Needs to be addressed
	if($newdStatus !~ /CI/ && $newpdStatus !~ /CI/){
	    print $rxn,"\tplantdefault\tstatus\t",$newpdStatus,"\n";

	    #If changes were made to the proton count in the equation, it's possible protons were added or removed
	    #In which case, the equation itself needs to be re-printed
	    if($newpdStatus =~ /HB/){
		print $rxn,"\tplantdefault\tequation\t",$pdRxn->equation(),"\n";
	    }
	}
    }

#    last;
    
}

__END__

#    print join("|", map { $_.":".$results->{$_} } keys %$results),"\t",$rObj->status(),"\n";
#    if(exists($results->{imbalancedAtoms})){
#	print "\t",join("|", map { $_.":".$results->{imbalancedAtoms}{$_}} keys %{$results->{imbalancedAtoms}}),"\n";
#    }
