#!/usr/bin/perl

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

#Handle filename
my $filestub=$file;
$filestub =~ s/\.xml$//;

my @compartments=@{$doc->getElementsByTagName('listOfCompartments')}[0]->getChildrenByTagName('*');
open(FH,"> ".$filestub."_Compartment_Table.txt");
foreach my $c (@compartments){
    print FH $c->getAttribute('id'),"\t",$c->getAttribute('name'),"\n";
}
close(FH);

my @compounds=@{$doc->getElementsByTagName('listOfSpecies')}[0]->getChildrenByTagName('*');
my ($id,$name,$cpt,$charge,$bC)=('','','',0,0);

open(FH,"> ".$filestub."_Compound_Table.txt");
print FH "ID\tName\tCompartment\tCharge\tboundary Condition\tKEGG\tInChI\n";
foreach my $c (@compounds){
    $id=$c->getAttribute('id');
    $name=$c->getAttribute('name');
    $cpt=$c->getAttribute('compartment');
    $charge=$c->getAttribute('charge');
    $bC=$c->getAttribute('boundaryCondition');

    if(!defined($id)){$id="UNDEF"};
    if(!defined($name)){$name="UNDEF"};
    if(!defined($cpt)){$cpt="UNDEF"};
    if(!defined($charge)){$charge="UNDEF"};
    if(!defined($bC)){$bC="UNDEF"};

    my @refs = $c->getElementsByTagName('rdf:li');
    my @KEGG=();
    my $InChI="";
    foreach my $ref (@refs){
	my $attr=$ref->getAttribute('rdf:resource');
	if($attr =~ /kegg/){
	    my @temp=split(/:/,$attr);
	    $attr=$temp[$#temp];
	    push(@KEGG,$attr);
	}
	if($attr =~ /inchi/){
	    my @temp=split(/\//,$attr);
	    $attr=join("/",@temp[4..$#temp]);
	    $InChI=$attr;
	}
    }
    print FH $id,"\t",$name,"\t",$cpt,"\t",$charge,"\t",$bC,"\t",join("|",@KEGG),"\t",$InChI,"\n";
}
close(FH);

my @notes=(); #CM
my @anns=(); #LN
my @reactions=@{$doc->getElementsByTagName('listOfReactions')}[0]->getChildrenByTagName('*');
my ($rev,$stoich,$cpd,$reacts,$prods)=(0,0,"","","");
my ($gene,$EC,$path)=("","","");
my @contents=();
my %genes=();
my %ECs=();
my %Paths=();

open(FH,"> ".$filestub."_Reaction_Table.txt");
print FH "ID\tName\tReversible\tReactants\tProducts\tE.C.\tPathway\n";

foreach my $c (@reactions){
    $id=$c->getAttribute('id');
    $name=$c->getAttribute('name');
    $rev=$c->getAttribute('reversible');    
    
    if(!defined($id)){$id="UNDEF"};
    if(!defined($name)){$name="UNDEF"};
    if(!defined($rev)){$rev="UNDEF"};
   
    @notes=$c->getChildrenByTagName('notes');
    if(scalar(@notes)>0){
	@notes=$notes[0]->getChildrenByTagName('html:p');

	%genes=();
	%ECs=();

	foreach my $n (@notes){
	    @contents=split(": ", $n->textContent);
	    if(substr($contents[0],0,4) eq "GENE"){
		my @genes=split(/or|and|,|snd/,$contents[1]);
		my $splice=undef;
		for(my $i=0;$i<scalar(@genes);$i++){
		    $genes[$i] =~ s/[\s;\(\)\[\]]//g;
		    if($genes[$i] eq ""){
			$splice=$i;
		    }
		}
		if(defined($splice)){
		    splice(@genes,$splice,1);
		}

		foreach my $g(@genes){
		    $genes{$g}=1;
		}
	    }elsif(substr($contents[0],0,4) eq "Prot"){
		$EC=$contents[1];
		$ECs{$EC}=1;
	    }else{
		print "Notes Warning: ",$n->textContent(),"\n";
	    }
	}
    }

    @anns=$c->getChildrenByTagName('annotation');
    if(scalar(@anns)>0){
	#four fields identified in LN models
	#EC, Gene Id, Gene Name, Pathway Name
        @anns=$anns[0]->getChildrenByTagName('html:p');

        %genes=();
        %ECs=();
	%Paths=();

        foreach my $n (@anns){
            @contents=split(": ", $n->textContent);
            if($contents[0] eq "Gene Id"){
                $gene=$contents[1];
                $genes{$gene}=1;
	    }elsif($contents[0] eq "Gene Name"){
                $gene=$contents[1];
                $genes{$gene}=1;
            }elsif($contents[0] eq "Pathway Name"){
                $path=$contents[1];
                $Paths{$gene}=1;
            }elsif($contents[0] eq "EC"){
                $EC=$contents[1];
                $ECs{$EC}=1;
	    }else{
		print "Annotation Warning: ",$n->textContent(),"\n";
	    }
	}
    }

    $reacts="";
    if(scalar(@{$c->getChildrenByTagName('listOfReactants')})>0){
	my @reactants_left=@{$c->getChildrenByTagName('listOfReactants')}[0]->getChildrenByTagName('*');
	foreach my $r (@reactants_left){
	    $stoich=$r->getAttribute('stoichiometry');
	    $cpd=$r->getAttribute('species');
	    
	    if(!defined($stoich)){$stoich="1"};
	    if(!defined($cpd)){$cpd="UNDEF"};
	    
	    $reacts.=$cpd."[".$stoich."];"
	}
    }
    
    $prods="";
    if(scalar(@{$c->getChildrenByTagName('listOfProducts')})>0){
	my @products_right=@{$c->getChildrenByTagName('listOfProducts')}[0]->getChildrenByTagName('*');
	foreach my $p (@products_right){
	    $stoich=$p->getAttribute('stoichiometry');
	    $cpd=$p->getAttribute('species');
	    if(!defined($stoich)){$stoich="1"};
	    if(!defined($cpd)){$cpd="UNDEF"};
	    
	    $prods.=$cpd."[".$stoich."];"
	}
    }

    ($gene,$EC,$path)=("","","");

    foreach my $ec (keys %ECs){
	if($ec eq ""){next;}
	$EC.=$ec.";";
    }
    foreach my $p (keys %Paths){
	if($p eq ""){next;}
	$path.=$p.";";
    }
    foreach my $g (keys %genes){
	if($g eq ""){next;}
	$gene.=$g.";";
    }

    print FH $id,"\t",$name,"\t",$rev,"\t",$reacts,"\t",$prods,"\t",$EC,"\t",$path,"\t",$gene,"\n";
 
}
close(FH);
