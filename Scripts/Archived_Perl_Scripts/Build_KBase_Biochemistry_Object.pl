#!/usr/bin/env perl
use Getopt::Long::Descriptive;
use warnings;
use strict;
use JSON;
my @temp=();
my $header=1;

my ($opt, $usage) = describe_options("%c %o <ID>",
	[ "compounds=s", "path to master compounds file", { default => "../Biochemistry/compounds.master.tsv" } ],
	[ "compartments=s", "path to master compartments file", { default => "../Biochemistry/compartments.master.tsv" } ],
	[ "reactions=s", "path to master reactions file", { default => "../Biochemistry/reactions.master.tsv" } ],
	[ "master=s", "path to output master biochemistry json file", { default => "../Biochemistry/kbase_biochemistry.master.json" } ],
	[ "help|h", "print usage message and exit" ]
    );

print($usage->text), exit if $opt->help;
my $id = $ARGV[0];
$id = 'master' if !$id;

# Create a new Biochemistry object (starting from scratch).
my $Biochem = {"id" => $id, "name"=> $id, description => "",
	       compartments => [],
	       compounds => [],
	       reactions => [],
	       reactionSets => [],
	       compoundSets => [],
	       cues => [],
	       compound_aliases => {},
	       reaction_aliases => {}};

#Not found in table: md5, unchargedFormula, cues
#Not found in object: is_core, is_obsolete, aliases, linked_compound
#Needs translation: structure -> structure_ref, abstract_compound -> abstractCompund_ref, comprised_of => comprisedOfCompound_refs

my %Compound_Headers=('is_cofactor' => 'isCofactor',
		      'charge' => 'defaultCharge',
		      'deltag' => 'deltaG',
		      'deltagerr' => 'deltaGErr',
		      'pka' => 'pkas',
		      'pkb' => 'pkbs');
my %Numeric_Headers = ('defaultCharge'=>1,'mass'=>1,'deltaG'=>1,'deltaGErr'=>1,'isCofactor'=>1);
my %Skipped_Headers = ('is_core'=>1,'is_obsolete'=>1,'aliases'=>1,'linked_compound'=>1);
my %Compounds=();
open(FH, "< ".$opt->compounds);
my @headers = map { exists($Compound_Headers{$_}) ? $Compound_Headers{$_} : $_ } split(/\t/,<FH>);
chomp $headers[$#headers];
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);

    my $cpdHash = {};
    for (my $i=0;$i<scalar(@temp);$i++) {
	next if exists($Skipped_Headers{$headers[$i]});

        if(exists($Numeric_Headers{$headers[$i]})){
	    $temp[$i]=0 if $temp[$i] eq 'null';
	    $temp[$i]+=0;
	}

	if($headers[$i] eq "structure"){

	}elsif($headers[$i] eq 'abstract_compound'){

	}elsif($headers[$i] eq 'comprised_of'){

	}elsif($headers[$i] =~ /^pk[ab]s$/){
	    next if $temp[$i] eq "null";
	    next if $temp[$i] =~ /^\s*$/;

	    my @pk = split(";", $temp[$i]);
	    my %pk=();

	    # Input argument is a hash keyed by dissocation constant values.
	    for (my $j=0; $j<scalar(@pk); $j++) {

		my @elements = split(":",$pk[$j]);
		$pk{$elements[0]}{$elements[1]}=1;
	    }
	    foreach my $number (keys %pk){
		$cpdHash->{$headers[$i]}{$number} = [ map { $_+=0 } sort { $a <=> $b } keys %{$pk{$number}} ];
	    }

	}else{
	    $cpdHash->{$headers[$i]}=$temp[$i];
	}
    }
    push(@{$Biochem->{compounds}},$cpdHash);
}
close(FH);
 
my $num_cpds = scalar(@{$Biochem->{compounds}});
print "Added ".$num_cpds." compounds to biochemistry object\n";

#Not found in table: md5, cues, defaultProtons
#Not found in object: code, equation, definition, aliases, linked_reaction, is_obsolete, compound_ids, pathways, is_transport, ec_numbers
#Needs translation: stoichiometry -> reagents, abstract_reaction -> abstractReaction_Ref

my %Reaction_Headers=('reversibility' => 'thermoReversibility',
		      'deltag' => 'deltaG',
		      'deltagerr' => 'deltaGErr');

%Numeric_Headers = ('deltaG'=>1,'deltaGErr'=>1);
%Skipped_Headers = ('code'=>1,'equation'=>1,'definition'=>1,'aliases'=>1,
		    'linked_reaction'=>1,'is_obsolete'=>1,'compound_ids'=>1,
		    'pathways'=>1,'is_transport'=>1,'ec_numbers'=>1);

my %Reactions=();
my %Compartments=();
open(FH, "< ".$opt->reactions);
@headers = map { exists($Reaction_Headers{$_}) ? $Reaction_Headers{$_} : $_ } split(/\t/,<FH>);
chomp $headers[$#headers];
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);

    my $rxnHash = {'defaultProtons'=>0};
    for (my $i=0;$i<scalar(@temp);$i++) {
	next if exists($Skipped_Headers{$headers[$i]});

        if(exists($Numeric_Headers{$headers[$i]})){
	    $temp[$i]=0 if $temp[$i] eq 'null';
	    $temp[$i]+=0;
	}

	if($headers[$i] eq 'abstract_reaction'){

	}elsif($headers[$i] eq 'stoichiometry'){
	    foreach my $rgt (split(/;/,$temp[$i])){
		my ($coefficient,$cpd,$cpt,$index,$name) = split(/:/,$rgt);
		$Compartments{$cpt}=1;
		my $rgtHash = {'compound_ref' => '~/compounds/id/'.$cpd,
			       'compartment_ref' => '~/compartments/id/'.$cpt,
			       'coefficient' => $coefficient+0,
			       'isCofactor' => 0};
		push(@{$rxnHash->{'reagents'}},$rgtHash);
	    }
	}else{
	    $rxnHash->{$headers[$i]}=$temp[$i];
	}
    }
    print "Error: ${temp[0]} doesn't have any reagents\n" if !exists($rxnHash->{'reagents'});
    next if !exists($rxnHash->{'reagents'});
    push(@{$Biochem->{reactions}},$rxnHash);
}

my $num_rxns = scalar(@{$Biochem->{reactions}});
print "Added ".$num_rxns." reactions to biochemistry object\n";

foreach my $cpt (sort keys %Compartments){
    my $cptHash = {'name'=>'generic','id'=>$cpt,'hierarchy'=>0};
    push(@{$Biochem->{compartments}},$cptHash);
}

my $num_cpts = scalar(@{$Biochem->{compartments}});
print "Added ".$num_cpts." compartments to biochemistry object\n";

open(OUT, "> ".$opt->master);
print OUT to_json($Biochem, { pretty => 1, ascii =>1 });
close(OUT);
