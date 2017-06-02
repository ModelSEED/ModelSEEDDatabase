#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my @temp=();
my $header = 1;

#typedef structure {
#    templatecompartment_id id;
#    string name;
#    list<string> aliases;
#    int hierarchy;
#    float pH;
#} TemplateCompartment;

open(CPT, "< ../Biochemistry/compartments.master.tsv");
my @Cpt_Headers = split(/\t/,<CPT>,-1);
chomp($Cpt_Headers[$#Cpt_Headers]);
my @Cpts=();
my %Cpts=();
while(<CPT>){
    chomp;
    @temp=split(/\t/,$_,-1);

    my $Cpt = {};
    for(my $i=0;$i<scalar(@Cpt_Headers);$i++){
	$Cpt->{$Cpt_Headers[$i]}=$temp[$i];
    }
    $Cpt->{aliases}=[];
    $Cpt->{pH}=7.0;

    push(@Cpts,$Cpt);
    $Cpts{$Cpt->{id}}=$Cpt;
}
close(CPT);

#typedef structure {
#    templatecompound_id id;
#    compound_ref compound_ref;
#    string name;
#    string abbreviation;
#    string md5;
#    bool isCofactor;
#    list<string> aliases;
#    float defaultCharge;
#    float mass;
#    float deltaG;
#    float deltaGErr;
#    string formula;
#} TemplateCompound;

open(CPD, "< ../Biochemistry/compounds.master.tsv");
my @Cpd_Headers = split(/\t/,<CPD>,-1);
chomp($Cpd_Headers[$#Cpd_Headers]);
my %Transform = ('is_cofactor' => 'isCofactor', 'charge' => 'defaultCharge', 'deltag' => 'deltaG', 'deltagerr' => 'deltaGErr');
my @Cpds=();
my %Cpds=();
while(<CPD>){
    chomp;
    @temp=split(/\t/,$_,-1);

    my $Cpd = {};
    for(my $i=0;$i<scalar(@Cpd_Headers);$i++){
	my $header = $Cpd_Headers[$i];
	next if $header =~ /^(source|structure|is_core|is_obsolete|linked_compound|pka|pkb|abstract_compound|comprised_of)$/;
	$header = exists($Transform{$header}) ? $Transform{$header} : $header;
	if($header eq "aliases"){
	    $Cpd->{$header}=[split(/\|/,$temp[$i])];
	}else{
	    $Cpd->{$header}=$temp[$i];
	}
	$Cpd->{compound_ref}="";
    }

    push(@Cpds,$Cpd);
    $Cpds{$Cpd->{id}}=$Cpd;
}
close(CPD);

open(RXN, "< ../Biochemistry/reactions.master.tsv");
my @Rxn_Headers = split(/\t/,<RXN>,-1);
chomp($Rxn_Headers[$#Rxn_Headers]);
my @Rxns=();
my %Rxns=();
while(<RXN>){
    chomp;
    @temp=split(/\t/,$_,-1);

    my $Rxn = {};
    for(my $i=0;$i<scalar(@Rxn_Headers);$i++){
	my $header = $Rxn_Headers[$i];
	$Rxn->{$header}=$temp[$i];
    }

    push(@Rxns,$Rxn);
    $Rxns{$Rxn->{id}}=$Rxn;
}
close(RXN);

#typedef structure {
#    modeltemplate_id id;
#    string name;
#    string modelType;
#    string domain;
#    Biochemistry_ref biochemistry_ref;

#    list<TemplateRole> roles;
#    list<TemplateComplex> complexes;
#    list<TemplateCompound> compounds;
#    list<TemplateCompCompound> compcompounds;
#    list<TemplateCompartment> compartments;
#    list<TemplateReaction> reactions;
#    list<TemplateBiomass> biomasses;
#    list<TemplatePathway> pathways;
#} ModelTemplate;

open(DATA, "< ../ModelTemplates/ModelTemplate_Data.txt");
my @Data_Headers = split(/\t/,<DATA>,-1);
chomp($Data_Headers[$#Data_Headers]);
while(<DATA>){
    chomp;
    @temp=split(/\t/,$_,-1);

    my $MT_Dir = "../ModelTemplates/".$temp[0]."_PMS/";
    next unless -d $MT_Dir;
    my $Map_Dir = "../Mappings/".$temp[6]."/";
    next unless -d $Map_Dir;

    my $ModelTemplate = {};
    for(my $i=0;$i<scalar(@Data_Headers);$i++){
	next if $Data_Headers[$i] =~ /^(ws_id|mapping_ref|mapping_id)$/;
	$ModelTemplate->{$Data_Headers[$i]}=$temp[$i];
    }

#    typedef structure {
#	templaterole_id id;
#	string name;
#	string source;
#	string type;
#	list<string> aliases;
#	list<feature_id> features;
#    } TemplateRole;

    open(ROLES, "< ".$Map_Dir."Mapping_Roles.txt");
    my @Role_Headers = split(/\t/,<ROLES>,-1);
    chomp($Role_Headers[$#Role_Headers]);
    my @Roles=();
    my %Roles=();
    while(<ROLES>){
	chomp;
	@temp=split(/\t/,$_,-1);

	my $Role = {};
	for(my $i=0;$i<scalar(@Role_Headers);$i++){
	    if($Role_Headers[$i] eq 'seedfeature'){
		$Role->{'features'} = [split(/;/,$temp[$i])];
	    }else{
		$Role->{$Role_Headers[$i]}=$temp[$i];
	    }
	}
	$Role->{source}="";
	$Role->{type}="";
	$Role->{aliases}=[];
	
	push(@Roles,$Role);
	$Roles{$Role->{id}}=$Role;
    }
    close(ROLES);

#    typedef structure {
#	templaterole_ref templaterole_ref;
#	int optionalRole;
#	string type;
#	int triggering;
#	float confidence;
#    } TemplateComplexRole;

#    typedef structure {
#	templatecomplex id;
#	string name;
#	string reference;
#	string source;
#	float confidence;
#	list<TemplateComplexRole> complexroles;
#    } TemplateComplex;

    open(CPXS, "< ".$Map_Dir."Mapping_Complexes.txt");
    my @Cpx_Headers = split(/\t/,<CPXS>,-1);
    chomp($Cpx_Headers[$#Cpx_Headers]);
    my @Cpxs=();
    my %Cpxs=();
    while(<CPXS>){
	chomp;
	@temp=split(/\t/,$_,-1);

	my $Cpx = {};
	for(my $i=0;$i<scalar(@Cpx_Headers);$i++){
	    if($Cpx_Headers[$i] eq 'complexroles'){
		foreach my $cpxrole (split(/\|/,$temp[$i])){
		    my $CpxRole = {};
		    foreach my $col (split(/\;/,$cpxrole)){
			my ($key,$value) = split(/:/,$col);
			if($key eq 'role_ref'){
			    $value =~ s/^\~\/roles\/id/~\/templateroles/;
			    $CpxRole->{templaterole_ref}=$value;
			}
			$CpxRole->{$key}=$value;
		    }
		    $CpxRole->{confidence}=1.0;
		    push(@{$Cpx->{complexroles}},$CpxRole);
		}
	    }else{
		$Cpx->{$Cpx_Headers[$i]}=$temp[$i];
	    }
	}
	$Cpx->{reference}="";
	$Cpx->{source}="";
	$Cpx->{confidence}=1.0;

	push(@Cpxs,$Cpx);
	$Cpxs{$Cpx->{id}}=$Cpx;
    }
    close(CPXS);

#    typedef structure {
#	templatereaction_id id;
#	reaction_ref reaction_ref;
#	string name;
#	string type;
#	string reference;
#	string direction;
#	string GapfillDirection;
#	float maxforflux;
#	float maxrevflux;
#	templatecompartment_ref templatecompartment_ref;
#	float base_cost;
#	float forward_penalty;
#	float reverse_penalty;
#	list<TemplateReactionReagent> templateReactionReagents;
#	list<templatecomplex_ref> templatecomplex_refs;
#    } TemplateReaction;

    open(TMPLRXN, "< ".$MT_Dir."ModelTemplate_Reactions.txt");
    my @TmplRxn_Headers = split(/\t/,<TMPLRXN>,-1);
    chomp($TmplRxn_Headers[$#TmplRxn_Headers]);
    my @TmplRxns=();
    my %TmplRxns=();
    while(<TMPLRXN>){
	chomp;
	@temp = split(/\t/,$_,-1);
	my $TmplRxn = {};
	my $Rxn = {};
	for(my $i=0;$i<scalar(@TmplRxn_Headers);$i++){
	    if($TmplRxn_Headers[$i] eq "reaction_ref"){
		my $rxn = substr($temp[$i],-8);
		$Rxn = $Rxns{$rxn};
	    }
	    
	    #Modify old reference
	    if($TmplRxn_Headers[$i] =~ /refs?$/){
		$temp[$i] =~ s/\/id\//\//g;
		$temp[$i] =~ s/\d+\/\d+\/\d+\//~\//g;
	    }

	    if($TmplRxn_Headers[$i] eq "compartment_ref"){
		$TmplRxn->{templatecompartment_ref}=$temp[$i];
	    }elsif($TmplRxn_Headers[$i] eq "complex_refs"){
		$TmplRxn->{templatecomplex_refs}=[split(/\|/,$temp[$i])];
	    }else{
		$TmplRxn->{$TmplRxn_Headers[$i]}=$temp[$i];
	    }
	}

#	typedef structure {
#	    templatecompcompound_id id;
#	    templatecompound_ref templatecompound_ref;
#	    float charge;
#	    float maxuptake;
#	    string formula;
#	    templatecompartment_ref templatecompartment_ref;
#	} TemplateCompCompound;

#	typedef structure {
#	    templatecompcompound_ref templatecompcompound_ref;
#	    float coefficient;
#	} TemplateReactionReagent;

	#Process equation
	my @TmplRxnRgts = ();
	foreach my $rgt (split(/;/,$Rxn->{stoichiometry})){
	    my ($coeff,$cpd,$cpt,$idx,$name) = split(/:/,$rgt);
	    my $TmplCptCpd = { id => '', 
			       templatecompound_ref => '', 
			       charge => $Cpds{$cpd}{defaultCharge},
			       maxuptake => 1000.0, 
			       formula => $Cpds{$cpd}{formula}, 
			       templatecompartment_ref => '' };

	    my $TmplRxnRgt = { templatecompcompound_ref => '',
			       coefficient => $coeff };

	    push(@TmplRxnRgts, $TmplRxnRgt);
	}
	$TmplRxn->{templateReactionReagents}=\@TmplRxnRgts;
	$TmplRxn->{maxforflux}=1000.0;
	$TmplRxn->{maxrevflux}=-1000.0;
	$TmplRxn->{reference}="~/templatereactions/".$TmplRxn->{id};

	$TmplRxns{$TmplRxn->{id}}=$TmplRxn;
	push(@TmplRxns,$TmplRxn);
    }
    close(TMPLRXN);

#    typedef structure {
#	templatebiomass_id id;
#	string name;
#	string type;
#	float other;
#	float dna;
#	float rna;
#	float protein;
#	float lipid;
#	float cellwall;
#	float cofactor;
#	float energy;
#	list<TemplateBiomassComponent> templateBiomassComponents;
#    } TemplateBiomass;

    open(TMPLBIO, "< ".$MT_Dir."ModelTemplate_Biomasses.txt");
    my @TmplBio_Headers = split(/\t/,<TMPLBIO>,-1);
    chomp($TmplBio_Headers[$#TmplBio_Headers]);
    my @TmplBios = ();
    while(<TMPLBIO>){
	chomp;
	@temp = split(/\t/,$_,-1);
	my $TmplBio = {};
	for(my $i=0;$i<scalar(@TmplBio_Headers);$i++){
	    $TmplBio->{$TmplBio_Headers[$i]}=$temp[$i];
	}

#	typedef structure {
#	    string class;
#	    templatecompcompound_ref templatecompcompound_ref;
#	    string coefficientType;
#	    float coefficient;
#	    list<templatecompcompound_ref> linked_compound_refs;
#	    list<float> link_coefficients;
#	} TemplateBiomassComponent;

	open(TMPLBIOCPN, "< ".$MT_Dir."ModelTemplate_Biomass_Components.txt");
	my @TmplBioCpn_Headers = split(/\t/,<TMPLBIOCPN>,-1);
	chomp($TmplBioCpn_Headers[$#TmplBioCpn_Headers]);
	my @TmplBioCpns=();
	while(<TMPLBIOCPN>){
	    chomp;
	    @temp = split(/\t/,$_,-1);

	    #Only biomass components whose id matches
	    my $TmplBio_ID = $TmplBio->{id};
	    next unless $temp[0] =~ /^${TmplBio_ID}/;

	    my $TmplBioCpn = {};
	    for(my $i=0;$i<scalar(@TmplBioCpn_Headers);$i++){
		next if $TmplBioCpn_Headers[$i] =~ /^(id|compound_ref|compartment_ref)$/;

		if($TmplBioCpn_Headers[$i] eq "linked_compound_refs" || $TmplBioCpn_Headers[$i] eq "link_coefficients"){
		    $TmplBioCpn->{$TmplBioCpn_Headers[$i]}=[split(/\|/,$temp[$i])];
		}else{
		    $TmplBioCpn->{$TmplBioCpn_Headers[$i]}=$temp[$i];
		}
	    }
	    $TmplBioCpn->{templatecompcompound_ref}="";

	    push(@TmplBioCpns,$TmplBioCpn);
	}
	$TmplBio->{templateBiomassComponents}=\@TmplBioCpns;
	push(@TmplBios,$TmplBio);
    }
    close(TMPLBIO);

    $ModelTemplate->{roles}=\@Roles;
    $ModelTemplate->{complexes}=\@Cpxs;

    #Need to remove compounds/compartments not found in reactions
    $ModelTemplate->{compounds}=\@Cpds;
    $ModelTemplate->{compartments}=\@Cpts;

    $ModelTemplate->{compcompounds}=[];
    $ModelTemplate->{reactions}=\@TmplRxns;
    $ModelTemplate->{biomasses}=\@TmplBios;
    $ModelTemplate->{pathways}=[];

    open(OUT, "> ".$MT_Dir."ModelTemplate.JSON");
    print OUT to_json($ModelTemplate,{pretty=>1});
    close(OUT);

}
close(DATA);
