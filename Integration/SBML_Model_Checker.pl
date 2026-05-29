#!/usr/bin/perl
use warnings;
use strict;
use XML::LibXML; 

my $file;

my $numArgs = $#ARGV + 1;
if ($numArgs eq 0) {
    print "Provide an SBML file please\n";
    exit;
}
$file  = $ARGV[0];

my $parser=XML::LibXML->new;
my $doc=$parser->parse_file($file);

print "Checking elements and attributes in document\n";
my @elements=$doc->getElementsByTagName('*');
my %established_elements=('compartment'=>1,
			  'listOfCompartments'=>1,
			  'listOfProducts'=>1,
			  'listOfReactants'=>1,
			  'listOfReactions'=>1,
			  'listOfSpecies'=>1,
			  'model'=>1,
			  'reaction'=>1,
			  'species'=>1,
			  'speciesReference'=>1,
    );

my @attributes=();
my %established_attributes=('charge'=>1,
			    'compartment'=>1,
			    'id'=>1,
			    'name'=>1,
			    'reversible'=>1,
			    'species'=>1,
			    'stoichiometry'=>1,
			    'boundaryCondition'=>1,
    );

my %checked={};
foreach my $e (@elements){
    if( ! exists($established_elements{$e->nodeName()}) and ! exists($checked{'e_'.$e->nodeName()})){
	print "Element Warning: ",$e->nodeName(),"\n";
	$checked{'e_'.$e->nodeName()}=1;
    }
    @attributes=$e->getAttributes();
    foreach my $a (@attributes){
	if( ! exists($established_attributes{$a->nodeName()}) and ! exists($checked{'e_'.$e->nodeName().'a_'.$a->nodeName()})){
	    print "Attribute Warning: ",$e->nodeName,":",$a->nodeName(),"\n";
	    $checked{'e_'.$e->nodeName().'a_'.$a->nodeName()}=1;
	}
    }
}
print "Finished checking elements and attributes in document\n";

my @cpds=$doc->getElementsByTagName("species");
my @rxns=$doc->getElementsByTagName("reaction");
print(scalar(@cpds)," compounds and ",scalar(@rxns)," reactions\n");
