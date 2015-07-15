#!/usr/bin/env perl
use warnings;
use strict;
use Bio::KBase::ObjectAPI::KBaseBiochem::Biochemistry;
use Data::Dumper;

# Create a new Biochemistry object (starting from scratch).  Set the reference so it
# is not in the default UUID format.
my $biochem = Bio::KBase::ObjectAPI::KBaseBiochem::Biochemistry->new({"id" => "master_v1.0.0", "_reference" => "~"});

# Note that the order of adding objects is important.

# Get the master list of compounds from the source file.
print "Getting compounds from ../Biochemistry/compounds.master.tsv ...\n";
my %Compounds=();
my @temp=();
open(FH, "< ../Biochemistry/compounds.master.tsv");
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
my $num_cpds = keys %Compounds;
print "Added ".$num_cpds." compounds to biochemistry object\n";

# Get the master list of compartments from the source file.
print "Getting compartments from ../Biochemistry/compartments.master.tsv\n";
my %Compartments=();
open(FH, "< ../Biochemistry/compartments.default.tsv");
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
print "Getting reactions from ../Biochemistry/reactions.master.tsv\n";
my %Reactions=();
open(FH, "< ../Biochemistry/reactions.master.tsv");
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
foreach my $rxn (sort keys %Reactions) {
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
}
my $num_rxns = keys %Reactions;
print "Added ".$num_rxns." reactions to biochemistry object\n";

# Get aliases from other sources.

# Save to a json format file.
print "Saving biochemistry to ../Biochemistry/biochemistry.master.json\n";
my $output = $biochem->export( { "format" => "json" } );
open(OUT, "> ../Biochemistry/biochemistry.master.json");
print OUT $output;

exit(0);
