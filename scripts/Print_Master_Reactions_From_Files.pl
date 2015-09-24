#!/usr/bin/env perl
use warnings;
use strict;
use Bio::KBase::ObjectAPI::KBaseBiochem::Biochemistry;
use Data::Dumper;
use Getopt::Long::Descriptive;

my ($opt, $usage) = describe_options("%c %o ",
	[ "compounds=s", "path to master compounds file", { default => "../Biochemistry/compounds.master.tsv" } ],
	[ "compartments=s", "path to master compartments file", { default => "../Biochemistry/compartments.master.tsv" } ],
	[ "reactions=s", "path to master reactions file", { default => "../Biochemistry/reactions.master.tsv" } ],
	[ "mods=s", "path to reaction modifications file", { default => "../Biochemistry/reactions.master.mods" } ],
	[ "help|h", "print usage message and exit" ]
);
print($usage->text), exit if $opt->help;

my @temp=();

#######################################################
#Initialization
#######################################################

# Get the master list of compounds from the source file.
# Note that the master list has already had modifications applied.
my %Compounds=();
open(FH, "< ".$opt->compounds);
my @headers = split(/\t/,<FH>);
chomp $headers[$#headers];
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    for (my $i=0;$i<scalar(@temp);$i++) { # Add all of the fields
	$Compounds{$temp[0]}{$headers[$i]} = $temp[$i];
    }
    
    # Remove unknown values from hash so compound object defaults are used instead.
    if ($Compounds{$temp[0]}{"charge"} eq "null") {
	delete $Compounds{$temp[0]}{"charge"};
    }
    if ($Compounds{$temp[0]}{"mass"} eq "null") {
	delete $Compounds{$temp[0]}{"mass"};
    }
    if ($Compounds{$temp[0]}{"deltag"} eq "null") {
	delete $Compounds{$temp[0]}{"deltag"};
    }
    if ($Compounds{$temp[0]}{"deltagerr"} eq "null") {
	delete $Compounds{$temp[0]}{"deltagerr"};
    }
    
    # Name needs to be specified as an array.
    $Compounds{$temp[0]}{"names"} = [ $Compounds{$temp[0]}{"name"} ];
    delete $Compounds{$temp[0]}{"structure"}; # do I need to build a CompoundStructure object here?

}
close(FH);

# Get the master list of compartments from the source file.
my %Compartments=();
open(FH, "< ".$opt->compartments);
@headers = split(/\t/,<FH>);
chomp $headers[$#headers];
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    for (my $i=0;$i<scalar(@temp);$i++) { # Add all of the fields
	$Compartments{$temp[0]}{$headers[$i]} = $temp[$i];
    }
}
close(FH);

# Load up the required modifications for reactions.
open(FH, "< ".$opt->mods);
my %Rxn_Mods=();
while(<FH>){
    chomp;
    @temp = split(/\t/,$_,-1);

    if($temp[2] eq "replace"){
	$Rxn_Mods{$temp[0]}{$temp[2]} = [$temp[3], $temp[4]];
    }elsif($temp[2] eq "coefficient"){
	$Rxn_Mods{$temp[0]}{$temp[2]}{$temp[3]}=$temp[4];
    }elsif($temp[2] eq "add"){
	$Rxn_Mods{$temp[0]}{$temp[2]}{$temp[3]}=[$temp[4], $temp[5]];
    }else{
	$Rxn_Mods{$temp[0]}{$temp[2]}=$temp[3];
    }
}
close(FH);

my %Rxns_Codes=();
my %Codes_Rxns=();
my %Reactions=(); #Master list of reactions

# Process each of the source files.  Note that if a reaction is defined in multiple
# source files, the reaction from the last file processed is used.
foreach my $db ("default", "plantdefault") {
    
    # The methods for building equations in a Reaction object can only be used if
    # there is a Biochemistry object with all of the valid compounds.  So first
    # create a new Biochemistry object and add all of the compounds and compartments
    # before processing reactions.

    # Create a new Biochemistry object (starting from scratch).  Set the reference so it
    # is not in the default UUID format.
    my $biochem = Bio::KBase::ObjectAPI::KBaseBiochem::Biochemistry->new({"id" => $db."_throw_away", "_reference" => "~"});

    # Add each compartment to the biochemistry object.
    foreach my $cpt (sort keys %Compartments) {
	# Add the compartment.
	$biochem->addCompartmentFromHash($Compartments{$cpt});
    }

    # Add each compartment to the biochemistry object.
    foreach my $cpd (sort keys %Compounds) {
	# Add the compound
	$biochem->addCompoundFromHash($Compounds{$cpd});
    }

    open(FH, "< ../Biochemistry/reactions.".$db.".tsv");
    @headers = split(/\t/,<FH>);
    chomp $headers[$#headers];
    while(<FH>){
	chomp;
	@temp=split(/\t/,$_,-1);
	
	# Assume that the id is in field 0.
	my $rxnid = $temp[0];

	#skip removed reactions
	next if exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{remove});
	
	# Skip the reaction when the modification is to give priority to different database.
	if (exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{priority}) && $Rxn_Mods{$rxnid}{priority} ne $db) {
	    next;
	}
	
	#skip the reaction if already retrieved from one database
	#this means that unless otherwise stated in the priority
	#we default to reactions from the default database
	next if exists($Reactions{$rxnid});
	
	my %rxnhash=(id=>$rxnid);
	for (my $i=1;$i<scalar(@temp);$i++) {
	    if (exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{$headers[$i]})) {
		# Apply a modification to a field.
		$temp[$i] = $Rxn_Mods{$rxnid}{$headers[$i]};
	    }
	    $rxnhash{$headers[$i]} = $temp[$i]; # Copy for adding to biochemistry object
	}
	
	# Remove unknown values so reaction object defaults are used instead.
	if ($rxnhash{"deltag"} eq "null") {
	    delete $rxnhash{"deltag"};
	}
	if ($rxnhash{"deltagerr"} eq "null") {
	    delete $rxnhash{"deltagerr"};
	}
	
	# Name needs to be specified as an array.
	$rxnhash{"names"} = [ $rxnhash{"name"} ];
	
	#update any modifications
	if (exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{replace})) {
	    my ($old_cpd,$new_cpd) = @{$Rxn_Mods{$rxnid}{replace}};

	    # Make sure replacement compound exists.
	    my $repCpd = $biochem->getObject("compounds", $old_cpd);
	    if (!$repCpd) {
		print "Replacement compound ".$old_cpd." does not exist\n";
		next;
	    }

	    $rxnhash{equation} =~ s/${old_cpd}/${new_cpd}/g;
	    $rxnhash{code} =~ s/${old_cpd}/${new_cpd}/g;
	    $rxnhash{stoichiometry} =~ s/${old_cpd}/${new_cpd}/g;

	    my $old_name = $Compounds{$old_cpd}{name};
	    my $new_name = $Compounds{$new_cpd}{name};

	    $rxnhash{definition} =~ s/${old_name}/${new_name}/g;
	}

	if (exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{coefficient})) {
	    foreach my $cpd (keys %{$Rxn_Mods{$rxnid}{coefficient}}){
		my $stoich = $Rxn_Mods{$rxnid}{coefficient}{$cpd};

		$rxnhash{equation} =~ s/\(\d+\) ${cpd}/\(${stoich}\) ${cpd}/;
		$rxnhash{code} =~ s/\(\d+\) ${cpd}/\(${stoich}\) ${cpd}/;

		$rxnhash{stoichiometry} =~ s/\d+:${cpd}/${stoich}:${cpd}/;

		my $name = $Compounds{$cpd}{name};
		$rxnhash{definition} =~ s/\(\d+\) ${name}/\(${stoich}\) ${name}/g;

	    }
	}

	if (exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{add})) {
	    foreach my $cpd (keys %{$Rxn_Mods{$rxnid}{add}}){
		my ($stoich,$cmpt) = @{$Rxn_Mods{$rxnid}{add}{$cpd}};

		if($stoich < 0){
		    #add compound to left side of equation
		    my $neutral_stoich = abs($stoich);
		    $rxnhash{equation} = "($neutral_stoich) ${cpd}[${cmpt}] ".$rxnhash{equation};
		    if($cpd ne 'cpd00067'){
			$rxnhash{equation} = "($neutral_stoich) ${cpd}[${cmpt}] ".$rxnhash{code};
		    }

		    my $name = $Compounds{$cpd}{name};
		    $rxnhash{stoichiometry} = "$neutral_stoich:$cpd:$cmpt:0:$name;".$rxnhash{stoichiometry};
		    $rxnhash{definition} = "($neutral_stoich) ${name}[${cmpt}] ".$rxnhash{definition};
		}elsif($stoich > 0){
		    #add compound to right side of equation
		    $rxnhash{equation} = $rxnhash{equation}." ($stoich) ${cpd}[${cmpt}]";
		    if($cpd ne 'cpd00067'){
			$rxnhash{equation} = $rxnhash{code}." ($stoich) ${cpd}[${cmpt}]";
		    }

		    my $name = $Compounds{$cpd}{name};
		    $rxnhash{stoichiometry} = $rxnhash{stoichiometry}.";$stoich:$cpd:$cmpt:0:$name";
		    $rxnhash{definition} = $rxnhash{definition}." ($stoich) ${name}[${cmpt}]";
		}
	    }
	}

	$rxnhash{findmatch}=0;
	my $rxn = $biochem->addReactionFromHash(\%rxnhash);
	$rxn->status($rxnhash{status});
	$rxn->checkReactionMassChargeBalance({rebalanceProtons=>1,rebalanceWater=>0,saveStatus=>1});
	
	# Generate a md5 of the equation for merging later.
	my $code = $rxn->genEquationCode();
	if (!exists($Codes_Rxns{$code})) { # If not unique, generate md5 with reversed equation.
	    $code = $rxn->revGenEquationCode();
	}
	
	# Map reaction IDs to md5.
	$Codes_Rxns{$code}{$rxnid}=1;
	$Rxns_Codes{$rxnid}{$code}=1;
	
	$Reactions{$rxnid}{equation} = $rxn->genEquation();
	$Reactions{$rxnid}{code} = $rxn->genCode();
	$Reactions{$rxnid}{definition} = $rxn->genDefinition();
	$Reactions{$rxnid}{stoichiometry} = $rxn->genStoichiometry();
	$Reactions{$rxnid}{status} = $rxn->status();
	
    }
    close(FH);
}

foreach my $db ("default", "plantdefault") {
    open(FH, "< ../Biochemistry/reactions.".$db.".tsv");
    @headers = split(/\t/,<FH>);
    chomp $headers[$#headers];
    while(<FH>){
	chomp;
	@temp=split(/\t/,$_,-1);
	
	# Assume that the id is in field 0.
	my $rxnid = $temp[0];
	
	#skip removed reactions
	next if exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{remove});

	# Skip the reaction when the modification is to give priority to different database.
	if (exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{priority}) && $Rxn_Mods{$rxnid}{priority} ne $db) {
	    next;
	}
	
	for (my $i=1;$i<scalar(@temp);$i++) {
	    if (exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{$headers[$i]})) {
		# Apply a modification to a field.
		$temp[$i] = $Rxn_Mods{$rxnid}{$headers[$i]};
	    }

	    if($headers[$i] !~ /status|equation|definition|code|stoichiometry/){
		$Reactions{$rxnid}{$headers[$i]} = $temp[$i];
	    }

	}
    }
}

#Print it all out
#avoiding re-visiting same reaction if already merged with another
my %Touched_Rxns=();
push(@headers,"is_obsolete");
push(@headers,"linked_reaction");
open(OUT, "> ".$opt->reactions);
print OUT join("\t",@headers),"\n",;
foreach my $rxn (sort keys %Rxns_Codes){

	#Obsolete if already 'touched' through merging
	$Reactions{$rxn}{is_obsolete} = (exists($Touched_Rxns{$rxn}) ? "1" : "0");
#	print Dumper($Reactions{$rxn});

	#Sorted so priority is given to lowest reaction identifier
	my %Merged_Rxns=();
	foreach my $code (sort keys %{$Rxns_Codes{$rxn}}){
		$Touched_Rxns{$rxn}{$code}=1;

		foreach my $merged_rxn (grep { $_ ne $rxn } sort keys %{$Codes_Rxns{$code}}) {
			$Merged_Rxns{$merged_rxn}=1;
			$Touched_Rxns{$merged_rxn}{$code}=1;
		}
	}

	$Reactions{$rxn}{linked_reaction} = (scalar(keys %Merged_Rxns)>0 ? join(";",sort keys %Merged_Rxns) : "null");

    #rebuild stoichiometry

#    print join("\n", map { $rxn.":".$_.":".$Rxns{$rxn}{$_} } grep { !defined($Rxns{$rxn}{$_}) } grep { $_ ne "id" } @headers),"\n";

	print OUT $rxn."\t".join("\t", map { $Reactions{$rxn}{$_} } grep { $_ ne "id" } @headers),"\n";
}
close(OUT);
