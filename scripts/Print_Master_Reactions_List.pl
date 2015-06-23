#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header=1;

#######################################################
#Initialization
#######################################################

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

#Load up the name modifications
#These are needed for the generation of stoichiometry field
open(FH, "< ../Biochemistry/compounds.master.mods");
my %Cpd_Mods=();
while(<FH>){
    chomp;
    @temp = split(/\t/,$_,-1);
    next unless $temp[2] eq "name";
    $Cpd_Mods{$temp[0]}{$temp[2]}=$temp[3];
}
close(FH);

#Load up the required modifications
open(FH, "< ../Biochemistry/reactions.master.mods");
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
my %Rxns=();

foreach my $db ("default","plantdefault"){
    my $bioObj = $FBAImpl->_get_msobject("Biochemistry","kbase",$db);

    #Modify compounds first
    foreach my $cpd (sort keys %Cpd_Mods){
	my $obj = $bioObj->getObject("compounds",$cpd);
	next if !$obj;

	foreach my $attr (keys %{$Cpd_Mods{$cpd}}){
	    $obj->$attr($Cpd_Mods{$cpd}{$attr});
	}
    }

    my @reactions = sort { $a->id() cmp $b->id() } grep { $_->id() ne "rxn00000" && $_->id() !~ /^CoA-2-methylpropanoylating/ } @{$bioObj->reactions()};
    foreach my $rxn (@reactions){

	if(exists($Rxn_Mods{$rxn->id()}) && exists($Rxn_Mods{$rxn->id()}{replace})){
	    my $repCpd = $bioObj->getObject("compounds", $Rxn_Mods{$rxn->id()}{replace}[1]);
	    next if !$repCpd;

	    foreach my $rgt (@{$rxn->reagents()}){
		if($rgt->compound()->id() eq $Rxn_Mods{$rxn->id()}{replace}[0]){
		    $rgt->compound($repCpd);
		}
	    }
	}

	next if exists($Rxn_Mods{$rxn->id()}) && exists($Rxn_Mods{$rxn->id()}{priority}) && $Rxn_Mods{$rxn->id()}{priority} ne $db;
	next if exists($Rxns{$rxn->id()});

	my $code = $rxn->genEquationCode();
	if(!exists($Codes_Rxns{$code})){
	    $code = $rxn->revGenEquationCode();
	}
	
	$Codes_Rxns{$code}{$rxn->id()}=1;
	$Rxns_Codes{$rxn->id()}{$code}=1;
	$Rxns{$rxn->id()}{equation}=$rxn->genEquation();
	$Rxns{$rxn->id()}{code}=$rxn->genCode();
	$Rxns{$rxn->id()}{definition}=$rxn->genDefinition();
	$Rxns{$rxn->id()}{stoichiometry}=$rxn->genStoichiometry();
    }
}

#Start with original biochemistry
open(FH, "< ../Biochemistry/reactions.default.tsv");
my @headers = split(/\t/,<FH>);
chomp $headers[$#headers];
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);

    if(exists($Rxn_Mods{$temp[0]}) && exists($Rxn_Mods{$temp[0]}{priority})){
	next;
    }

    for(my $i=1;$i<scalar(@temp);$i++){

	if(exists($Rxn_Mods{$temp[0]}) && exists($Rxn_Mods{$temp[0]}{$headers[$i]})){
	    $temp[$i] = $Rxn_Mods{$temp[0]}{$headers[$i]};
	}elsif($headers[$i] =~ /equation|definition|code|stoichiometry/){
	    next;
	}

	$Rxns{$temp[0]}{$headers[$i]}=$temp[$i];

    }
}
close(FH);

#foreach my $rxn (sort keys %Rxns_Codes){
#    print join("\n", map { $rxn.":".$_.":".$Rxns{$rxn}{$_} } grep { !defined($Rxns{$rxn}{$_}) } grep { $_ ne "id" } @headers),"\n";
#    last;
#}

#Add new biochemistry
open(FH, "< ../Biochemistry/reactions.plantdefault.tsv");
$header=1;
while(<FH>){
    chomp;
    if($header){$header--;next;}
    @temp=split(/\t/,$_,-1);
    next if exists($Rxns{$temp[0]}) && !exists($Rxns{$temp[0]}{equation});

    for(my $i=1;$i<scalar(@temp);$i++){
	if(exists($Rxn_Mods{$temp[0]}) && exists($Rxn_Mods{$temp[0]}{$headers[$i]})){
	    $temp[$i] = $Rxn_Mods{$temp[0]}{$headers[$i]};
	}elsif($headers[$i]  =~ /equation|definition|code|stoichiometry/){
	    next;
	}

	$Rxns{$temp[0]}{$headers[$i]}=$temp[$i];

    }
}
close(FH);

#Print it all out
#avoiding re-visiting same reaction if already merged with another
my %Touched_Rxns=();
push(@headers,"is_obsolete");
push(@headers,"linked_reaction");
open(OUT, "> Master_Reaction_List.tsv");
print OUT join("\t",@headers),"\n",;
foreach my $rxn (sort keys %Rxns_Codes){

    #Obsolete if already 'touched' through merging
    $Rxns{$rxn}{is_obsolete} = (exists($Touched_Rxns{$rxn}) ? "1" : "0");

    #Sorted so priority is given to lowest reaction identifier
    my %Merged_Rxns=();
    foreach my $code (sort keys %{$Rxns_Codes{$rxn}}){
	$Touched_Rxns{$rxn}{$code}=1;

	foreach my $merged_rxn (grep { $_ ne $rxn } sort keys %{$Codes_Rxns{$code}}){
	    $Merged_Rxns{$merged_rxn}=1;
	    $Touched_Rxns{$merged_rxn}{$code}=1;
	}
    }

    $Rxns{$rxn}{linked_reaction} = (scalar(keys %Merged_Rxns)>0 ? join(";",sort keys %Merged_Rxns) : "null");

    #rebuild stoichiometry

#    print join("\n", map { $rxn.":".$_.":".$Rxns{$rxn}{$_} } grep { !defined($Rxns{$rxn}{$_}) } grep { $_ ne "id" } @headers),"\n";

    print OUT $rxn."\t".join("\t", map { $Rxns{$rxn}{$_} } grep { $_ ne "id" } @headers),"\n";
}
close(OUT);
