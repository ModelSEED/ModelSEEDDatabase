#!/usr/bin/env perl
use warnings;
use strict;
use Bio::KBase::ObjectAPI::KBaseBiochem::Biochemistry;
use Bio::P3::Workspace::WorkspaceClientExt;
use Data::Dumper;
use Getopt::Long::Descriptive;

my ($opt, $usage) = describe_options("%c %o <ID>",
	[ "compounds=s", "path to master compounds file", { default => "../Biochemistry/compounds.master.tsv" } ],
	[ "compartments=s", "path to master compartments file", { default => "../Biochemistry/compartments.master.tsv" } ],
	[ "reactions=s", "path to master reactions file", { default => "../Biochemistry/reactions.master.tsv" } ],
	[ "master=s", "path to output master biochemistry json file", { default => "../Biochemistry/biochemistry.master.json" } ],
	[ "wsobject=s", "path to workspace object", { default => undef } ],
	[ "help|h", "print usage message and exit" ]
);
print($usage->text), exit if $opt->help;
my $id = $ARGV[0];

# Create a new Biochemistry object (starting from scratch).  Set the reference so it
# is not in the default UUID format.
my $biochem = Bio::KBase::ObjectAPI::KBaseBiochem::Biochemistry->new({"id" => $id, "_reference" => "~"});

# Note that the order of adding objects is important.

# Get the master list of compounds from the source file.
print "Getting compounds from ".$opt->compounds." ...\n";
my %Compounds=();
my @temp=();
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
	
	# Convert pkas and pkbs from file format to input argument.
	if ($Compounds{$cpd}{"pka"} eq "null") {
		delete $Compounds{$cpd}{"pka"};
	} else {
		# Multiple pkas are separated by semicolon.
		my @pka = split(";", $Compounds{$cpd}{"pka"});
		# Input argument is a hash keyed by dissocation constant values.
		$Compounds{$cpd}{"pkas"} = {};
		for (my $i=0; $i<scalar(@pka); $i++) {
			# Each pka is the format "atoms:value".
			my @elements = split(":",$pka[$i]);
			$Compounds{$cpd}{"pkas"}{$elements[1]} = [ int($elements[0]) ]; # Assuming only two elements after split
		}
		delete $Compounds{$cpd}{"pka"};
	}
	if ($Compounds{$cpd}{"pkb"} eq "null") {
		delete $Compounds{$cpd}{"pkb"};
	} else {
		# pkbs are handled the same way as pkas.
		my @pkb = split(";", $Compounds{$cpd}{"pkb"});
		$Compounds{$cpd}{"pkbs"} = {};
		for (my $i=0; $i<scalar(@pkb); $i++) {
			my @elements = split(":",$pkb[$i]);
			$Compounds{$cpd}{"pkbs"}{$elements[1]} = [ int($elements[0]) ]; # Assuming only two elements after split
		}
		delete $Compounds{$cpd}{"pkb"};
	}
	
	delete $Compounds{$cpd}{"structure"}; # do I need to build a CompoundStructure object here?
	
	# Add the compound.
	$biochem->addCompoundFromHash($Compounds{$cpd});
}
my $num_cpds = keys %Compounds;
print "Added ".$num_cpds." compounds to biochemistry object\n";

# Get the master list of compartments from the source file.
print "Getting compartments from ".$opt->compartments." \n";
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
my $num_cpts = keys %Compartments;
print "Added ".$num_cpts." compartments to biochemistry object\n";

# Get the master list of reactions from the source file.
print "Getting reactions from ".$opt->reactions."\n";
my %Reactions=();
open(FH, "< ".$opt->reactions);
@headers = split(/\t/,<FH>);
chomp $headers[$#headers];
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    for (my $i=0;$i<scalar(@temp);$i++) { # Add all of the fields
		$Reactions{$temp[0]}{$headers[$i]} = $temp[$i];
    }
}
close(FH);

# Add each reaction to the biochemistry object.
my $num_rxns = 0;
foreach my $rxn (sort keys %Reactions) {
	# Skip obsolete reactions.
	if ($Reactions{$rxn}{"is_obsolete"} eq "1") {
		next;
	}

	# Remove unknown values so reaction object defaults are used instead.
	if ($Reactions{$rxn}{"deltag"} eq "null") {
		delete $Reactions{$rxn}{"deltag"};
	}
	if ($Reactions{$rxn}{"deltagerr"} eq "null") {
		delete $Reactions{$rxn}{"deltagerr"};
	}

	# Name needs to be specified as an array.
	$Reactions{$rxn}{"names"} = [ $Reactions{$rxn}{"name"} ];
	
	# Add the reaction.
	$biochem->addReactionFromHash($Reactions{$rxn});
	$num_rxns++;
}
print "Added ".$num_rxns." reactions to biochemistry object\n";

# Get aliases from other sources.

# Save to a json format file.
print "Saving biochemistry to ".$opt->master."\n";
my $data = $biochem->export( { "format" => "json" } );
open(OUT, "> ".$opt->master);
print OUT $data;

if (defined($opt->wsobject)) {
	print "Storing biochemistry in workspace object ".$opt->wsobject."\n";
	my $wsClient = Bio::P3::Workspace::WorkspaceClientExt->new();
	my $output = $wsClient->create( { "objects" => [ [ $opt->wsobject, 'biochemistry', {}, $data ] ] });
}
exit(0);
