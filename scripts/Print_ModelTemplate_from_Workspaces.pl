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
                                                       'mfatoolkitbin' => "/homes/chenry/kbase/MFAToolkit/bin/mfatoolkit",
                                                       'probanno-url' => "http://140.221.85.86:7073/",
                                                       'mssserver-url' => "http://bio-data-1.mcs.anl.gov/services/ms_fba",
                                                       'accounttype' => "kbase",
                                                       'workspace-url' => "http://kbase.us/services/ws",
                                                       'defaultJobState' => "queued",
                                                       'gaserver-url' => "http://kbase.us/services/genome_annotation",
                                                       'idserver-url' => "http://kbase.us/services/idserver"});
$FBAImpl->_setContext(undef,{auth=>$AToken});
my $WSClient = $FBAImpl->_workspaceServices();

my @TmplRxn_Headers = ('id', 'reaction_ref', 'compartment_ref', 'complex_refs', 'direction', 'GapfillDirection', 'type', 'base_cost', 'forward_penalty', 'reverse_penalty');
my @TmplBio_Headers = ('id', 'name', 'type', 'other', 'dna', 'rna', 'protein', 'lipid', 'cellwall', 'cofactor', 'energy');
my @TmplBioCmp_Headers = ('id', 'class', 'compound_ref', 'compartment_ref', 'coefficientType', 'coefficient', 'linked_compound_refs', 'link_coefficients');

my @Data_Headers = ('id', 'name', 'modelType', 'domain', 'mapping_ref', 'biochemistry_ref');
open(DATA, '> ../ModelTemplates/ModelTemplate_Data.txt');
print DATA "ws_id\t",join("\t", @Data_Headers),"\n";

my %Mappings = ();
foreach my $template (@{$WSClient->list_objects({workspaces=>["KBaseTemplateModels"],type=>"KBaseFBA.ModelTemplate"})}){
    print $template->[1],"\n";
    my $MT_Dir = '../ModelTemplates/'.$template->[1].'/';
    mkdir $MT_Dir if !-d $MT_Dir;

    my $templateObj = $FBAImpl->_get_msobject("ModelTemplate","KBaseTemplateModels",$template->[1]);

    #Retrieve mappings
    my $mappingObj = (@{$WSClient->get_objects([{ref=>$templateObj->mapping_ref()}])})[0];
    my $mapping_wsid = $mappingObj->{info}[1];

    $mappingObj = Bio::KBase::ObjectAPI::KBaseOntology::Mapping->new($mappingObj->{data});
    $Mappings{$mapping_wsid}=$mappingObj;

    #Print out top-level data
    print DATA $template->[1],"\t";
    foreach my $header (@Data_Headers){
	print DATA $templateObj->$header();
	print DATA $header eq $Data_Headers[$#Data_Headers] ? "\n" : "\t";
    }

    open(TMPLRXN, "> ".$MT_Dir."ModelTemplate_Reactions.txt");
    print TMPLRXN join("\t", @TmplRxn_Headers),"\n";
    foreach my $tmplrxn (@{$templateObj->templateReactions()}){
	foreach my $header (@TmplRxn_Headers){
	    if($header eq 'complex_refs'){
		print TMPLRXN join("|", sort @{$tmplrxn->$header()});
	    }else{
		print TMPLRXN defined($tmplrxn->$header()) ? $tmplrxn->$header() : "";
	    }

	    print TMPLRXN $header eq $TmplRxn_Headers[$#TmplRxn_Headers] ? "\n" : "\t";
	}
    }
    close(TMPLRXN);

    open(TMPLBIO, "> ".$MT_Dir."ModelTemplate_Biomasses.txt");
    print TMPLBIO join("\t", @TmplBio_Headers),"\n";
    open(TMPLBIOCMP, "> ".$MT_Dir."ModelTemplate_Biomass_Components.txt");
    print TMPLBIOCMP join("\t", @TmplBioCmp_Headers),"\n";
    foreach my $tmplbio (@{$templateObj->templateBiomasses()}){
	foreach my $header (@TmplBio_Headers){
	    print TMPLBIO defined($tmplbio->$header()) ? $tmplbio->$header() : "";
	    print TMPLBIO $header eq $TmplBio_Headers[$#TmplBio_Headers] ? "\n" : "\t";
	}
	foreach my $tmplbiocmp (@{$tmplbio->templateBiomassComponents()}){
	    foreach my $header (@TmplBioCmp_Headers){
		if($header eq 'linked_compound_refs' || $header eq 'link_coefficients'){
		    print TMPLBIOCMP join("|", sort @{$tmplbiocmp->$header()});
		}else{
		    print TMPLBIOCMP defined($tmplbiocmp->$header()) ? $tmplbiocmp->$header() : "";		    
		}
		print TMPLBIOCMP $header eq $TmplBioCmp_Headers[$#TmplBioCmp_Headers] ? "\n" : "\t";
	    }
	}
    }
    close(TMPLBIO);
    close(TMPLBIOCMP);
}
close(DATA);

my @Roles_Headers = ('id', 'name', 'seedfeature');
my @Subsystems_Headers = ('id', 'name', 'class', 'subclass', 'type', 'role_refs');
my @Complexes_Headers = ('id', 'name', 'complexroles');
@Data_Headers = ('id', 'name', 'roles', 'subsystems', 'complexes', 'role_aliases', 'complex_aliases', 'subsystem_aliases');
open(DATA, '> ../Mappings/Mapping_Data.txt');
print DATA "ws_id\t",join("\t", @Data_Headers),"\n";
foreach my $mapping (keys %Mappings){
    my $Map_Dir = '../Mappings/'.$mapping.'/';
    mkdir $Map_Dir if !-d $Map_Dir;

    #Print out top-level data
    print DATA $mapping,"\t";
    foreach my $header (@Data_Headers){
#	if($header eq 'roles' || $header eq 'subsystems' || $header eq 'complexes'){
#	    print DATA join("|", sort @{$Mappings{$mapping}->$header()});
#	}else{
	    print DATA $Mappings{$mapping}->$header();
#	}

	print DATA $header eq $Data_Headers[$#Data_Headers] ? "\n" : "\t";
    }

    open(ROLES, "> ".$Map_Dir."Mapping_Roles.txt");
    print ROLES join("\t", @Roles_Headers),"\n";
    foreach my $role (@{$Mappings{$mapping}->roles()}){
	foreach my $header (@Roles_Headers){
	    print ROLES defined($role->$header()) ? $role->$header() : "";
	    print ROLES $header eq $Roles_Headers[$#Roles_Headers] ? "\n" : "\t";
	}
    }
    close(ROLES);

    open(COMPLEXES, "> ".$Map_Dir."Mapping_Complexes.txt");
    print COMPLEXES join("\t", @Complexes_Headers),"\n";
    foreach my $complex (@{$Mappings{$mapping}->complexes()}){
	foreach my $header (@Complexes_Headers){
	    if($header eq 'complexroles'){
		print COMPLEXES join("|", map { $_->role_ref().";".$_->optionalRole().";".$_->type().";".$_->triggering()} @{$complex->complexroles()});
	    }else{
		print COMPLEXES defined($complex->$header()) ? $complex->$header() : "";
	    }
	    print COMPLEXES $header eq $Complexes_Headers[$#Complexes_Headers] ? "\n" : "\t";
	}
    }
    close(COMPLEXES);

    open(SUBSYSTEMS, "> ".$Map_Dir."Mapping_Subsystems.txt");
    print SUBSYSTEMS join("\t", @Subsystems_Headers),"\n";
    foreach my $subsystem (@{$Mappings{$mapping}->subsystems()}){
	foreach my $header (@Subsystems_Headers){
	    if($header eq 'role_refs'){
		print SUBSYSTEMS join("|", sort @{$subsystem->$header()});
	    }else{
		print SUBSYSTEMS defined($subsystem->$header()) ? $subsystem->$header() : "";
	    }
	    print SUBSYSTEMS $header eq $Subsystems_Headers[$#Subsystems_Headers] ? "\n" : "\t";
	}
    }
    close(SUBSYSTEMS);    
}
close(DATA);
