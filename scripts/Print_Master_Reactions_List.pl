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

# The methods for building equations in a Reaction object can only be used if
# there is a Biochemistry object with all of the valid compounds.  So first
# create a new Biochemistry object and add all of the compounds and compartments
# before processing reactions.

# Create a new Biochemistry object (starting from scratch).  Set the reference so it
# is not in the default UUID format.
my $biochem = Bio::KBase::ObjectAPI::KBaseBiochem::Biochemistry->new({"id" => "throw_away", "_reference" => "~"});

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
}
close(FH);

# Add each compound to the biochemistry object.
foreach my $cpd (sort keys %Compounds) { # grep { $_ ne "cpd00000" }
	# Remove unknown values from hash so compound object defaults are used instead.
	if ($Compounds{$cpd}{"charge"} eq "null") {
		delete $Compounds{$cpd}{"charge"};
	}
	if ($Compounds{$cpd}{"mass"} eq "null") {
		delete $Compounds{$cpd}{"mass"};
	}
	if ($Compounds{$cpd}{"deltag"} eq "null") {
		delete $Compounds{$cpd}{"deltag"};
	}
	if ($Compounds{$cpd}{"deltagerr"} eq "null") {
		delete $Compounds{$cpd}{"deltagerr"};
	}
	
	# Name needs to be specified as an array.
	$Compounds{$cpd}{"names"} = [ $Compounds{$cpd}{"name"} ];
	delete $Compounds{$cpd}{"structure"}; # do I need to build a CompoundStructure object here?
	
	# Add the compound.
	$biochem->addCompoundFromHash($Compounds{$cpd});
}

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

# Add each compartment to the biochemistry object.
foreach my $cpt (sort keys %Compartments) {
	# Add the compartment.
	$biochem->addCompartmentFromHash($Compartments{$cpt});
}

# Load up the required modifications for reactions.
open(FH, "< ".$opt->mods);
my %Rxn_Mods=();
while(<FH>){
	chomp;
	@temp = split(/\t/,$_,-1);
	$Rxn_Mods{$temp[0]}{$temp[2]}=$temp[3];
	if($temp[2] eq "replace"){
		$Rxn_Mods{$temp[0]}{$temp[2]} = [$temp[3], $temp[4]];
	}
}
close(FH);

my %Rxns_Codes=();
my %Codes_Rxns=();
my %Reactions=(); # Master list of reactions

# Process each of the source files.  Note that if a reaction is defined in multiple
# source files, the reaction from the last file processed is used.
foreach my $db ("default", "plantdefault") {
	open(FH, "< ../Biochemistry/reactions.".$db.".tsv");
	@headers = split(/\t/,<FH>);
	chomp $headers[$#headers];
	while(<FH>){
	    chomp;
	    @temp=split(/\t/,$_,-1);
	
		# Assume that the id is in field 0.
		my $rxnid = $temp[0];
		my %rxnhash=();
		# Skip the reaction when the modification is to give priority to different database.
		if (exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{priority}) && $Rxn_Mods{$rxnid}{priority} ne $db) {
			next;
		}
	
		for (my $i=1;$i<scalar(@temp);$i++) {
			if (exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{$headers[$i]})) {
				# Apply a modification to a field.
				$temp[$i] = $Rxn_Mods{$rxnid}{$headers[$i]};
			}
	
			# Add this field to the reaction in the master list.
			$Reactions{$rxnid}{$headers[$i]} = $temp[$i];
			$rxnhash{$rxnid}{$headers[$i]} = $temp[$i]; # Copy for adding to biochemistry object
		}
	
		# Remove unknown values so reaction object defaults are used instead.
		if ($rxnhash{$rxnid}{"deltag"} eq "null") {
			delete $rxnhash{$rxnid}{"deltag"};
		}
		if ($rxnhash{$rxnid}{"deltagerr"} eq "null") {
			delete $rxnhash{$rxnid}{"deltagerr"};
		}
	
		# Name needs to be specified as an array.
		$rxnhash{$rxnid}{"names"} = [ $rxnhash{$rxnid}{"name"} ];
		
		# Add the reaction.
		$rxnhash{$rxnid}{"id"} = $rxnid;
		my $rxn = $biochem->addReactionFromHash($rxnhash{$rxnid});
	
		# Apply modification to replace a compound.
		if (exists($Rxn_Mods{$rxnid}) && exists($Rxn_Mods{$rxnid}{replace})) {
			# Make sure replacement compound exists.
			my $repCpd = $biochem->getObject("compounds", $Rxn_Mods{$rxnid}{replace}[1]);
			if (!$repCpd) {
				print "Replacement compound ".$Compounds{$Rxn_Mods{$rxnid}{replace}[1]}." does not exist\n";
				next;
			}
	
			# Find the compound in the list of reagents and replace it.
			my $foorxn = $biochem->searchForReaction($rxnid);
			foreach my $rgt (@{$foorxn->reagents()}) {
				if ($rgt->compound()->id() eq $Rxn_Mods{$rxnid}{replace}[0]) {
					$rgt->compound($repCpd);
				}
			}
		}
	
		# Generate a md5 of the equation for merging later.
		my $code = $rxn->genEquationCode();
		if (!exists($Codes_Rxns{$code})) { # If not unique, generate md5 with reversed equation.
			$code = $rxn->revGenEquationCode();
		}
	
		# Map reaction IDs to md5.
		$Codes_Rxns{$code}{$rxnid}=1;
		$Rxns_Codes{$rxnid}{$code}=1;
		
		# Update fields that use compounds.
		$Reactions{$rxnid}{equation} = $rxn->genEquation();
		$Reactions{$rxnid}{code} = $rxn->genCode();
		$Reactions{$rxnid}{definition} = $rxn->genDefinition();
		$Reactions{$rxnid}{stoichiometry} = $rxn->genStoichiometry();
		
	}
	close(FH);
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
