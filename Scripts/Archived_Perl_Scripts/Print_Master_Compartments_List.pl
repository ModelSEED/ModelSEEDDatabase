#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long::Descriptive;

my ($opt, $usage) = describe_options("%c %o ",
	[ "compartments=s", "path to master compartments file", { default => "../Biochemistry/compartments.master.tsv" } ],
	[ "help|h", "print usage message and exit" ]
);
print($usage->text), exit if $opt->help;

my @temp=();
my @headers=();
my %Compartments=(); # Master list of compartments

# Process each of the source files.  Note that if a compartment is defined in multiple
# source files, the compartment from the last file processed is used.
foreach my $db ("default", "plantdefault") {
	open(FH, "< ../Biochemistry/compartments.".$db.".tsv");
	@headers = split(/\t/,<FH>);
	chomp $headers[$#headers];
	while(<FH>){
	    chomp;
	    @temp=split(/\t/,$_,-1);

		for (my $i=1;$i<scalar(@temp);$i++) {
			$Compartments{$temp[0]}{$headers[$i]}=$temp[$i]; # Assume that the id is in field 0
		}
	}
}

# Create the master compartments file.
open(OUT, "> ".$opt->compartments);
print OUT join("\t",@headers),"\n";
foreach my $cpt (sort keys %Compartments) {
	print OUT $cpt."\t".join("\t", map { $Compartments{$cpt}{$_} } grep { $_ ne "id" } @headers),"\n";
}
