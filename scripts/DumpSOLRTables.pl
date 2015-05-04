#!/usr/bin/env perl
use warnings;
use strict;

#######################################################
#Initialization
#######################################################

my $directory = $ARGV[0];

use Bio::KBase::fbaModelServices::ScriptHelpers qw( getToken );
my $AToken = getToken();

use Bio::KBase::fbaModelServices::Impl;
my $fba = Bio::KBase::fbaModelServices::Impl->new({'fbajobcache' => "/homes/seaver/Projects/KBase_Scripts/FBA_Scripts/JobsCache",
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
$fba->_setContext(undef,{auth=>$AToken});
my $bio = $fba->_get_msobject("Biochemistry","kbase","default");
my $pbio = $fba->_get_msobject("Biochemistry","kbase","plantdefault");
my $map = $fba->_get_msobject("Mapping","kbase","default-mapping");

my $rxn_alias_hash = $bio->reactionsByAlias();
my $p_rxn_alias_hash = $pbio->reactionsByAlias();
my $rxn_pathways = {};
open(my $fh, "<", "/Users/chenry/workspace/maindb/HopeScenarios.txt");
while (my $line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	my $patharray = [split(/:/,$array->[0])];
	pop(@{$patharray});
	pop(@{$patharray});
	if (defined($rxn_alias_hash->{KEGG}->{$array->[1]})) {
		foreach my $rxn (keys(%{$rxn_alias_hash->{KEGG}->{$array->[1]}})) {
			$rxn_pathways->{$rxn}->{KEGG}->{$patharray->[0]} = 1;
			$rxn_pathways->{$rxn}->{Scenario}->{join("/",@{$patharray})} = 1;
		}
	} elsif (defined($p_rxn_alias_hash->{KEGG}->{$array->[1]})) {
		foreach my $rxn (keys(%{$p_rxn_alias_hash->{KEGG}->{$array->[1]}})) {
			$rxn_pathways->{$rxn}->{KEGG}->{$patharray->[0]} = 1;
			$rxn_pathways->{$rxn}->{Scenario}->{join("/",@{$patharray})} = 1;
		}
	}
}
close($fh);

my $cpd_alias_hash = $bio->compoundsByAlias();
my $p_cpd_alias_hash = $pbio->compoundsByAlias();
my $cpd_structure = {};
open($fh, "<", "/Users/chenry/code/ModelSEEDDatabase/Structures/KEGG_Search_InChI.txt");
while (my $line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	if (defined($cpd_alias_hash->{KEGG}->{$array->[0]})) {
		my $idarray = [keys(%{$cpd_alias_hash->{KEGG}->{$array->[0]}})];
		my $cpdid = $idarray->[0];
		if (!defined($cpd_structure->{$cpdid})) {
			$cpd_structure->{$cpdid} = $array->[1];
		}
	} elsif (defined($p_cpd_alias_hash->{KEGG}->{$array->[0]})) {
		my $idarray = [keys(%{$p_cpd_alias_hash->{KEGG}->{$array->[0]}})];
		my $cpdid = $idarray->[0];
		if (!defined($cpd_structure->{$cpdid})) {
			$cpd_structure->{$cpdid} = $array->[1];
		}
	}
}
close($fh);
open($fh, "<", "/Users/chenry/code/ModelSEEDDatabase/Structures/MetaCyc_Search_InChI.txt");
while (my $line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	if (defined($cpd_alias_hash->{MetaCyc}->{$array->[0]})) {
		my $idarray = [keys(%{$cpd_alias_hash->{KEGG}->{$array->[0]}})];
		my $cpdid = $idarray->[0];
		if (!defined($cpd_structure->{$cpdid})) {
			$cpd_structure->{$cpdid} = $array->[1];
		}
	} elsif (defined($p_cpd_alias_hash->{MetaCyc}->{$array->[0]})) {
		my $idarray = [keys(%{$p_cpd_alias_hash->{MetaCyc}->{$array->[0]}})];
		my $cpdid = $idarray->[0];
		if (!defined($cpd_structure->{$cpdid})) {
			$cpd_structure->{$cpdid} = $array->[1];
		}
	}
}
close($fh);

#Retreiving templates
my $templates = [
	$fba->_get_msobject("ModelTemplate","KBaseTemplateModels","GramPosModelTemplate"),
	$fba->_get_msobject("ModelTemplate","KBaseTemplateModels","GramNegModelTemplate"),
	$fba->_get_msobject("ModelTemplate","KBaseTemplateModels","CoreModelTemplate"),
	$fba->_get_msobject("ModelTemplate","KBaseTemplateModels","PlantModelTemplate")
];
#Printing template table
my $templatlist = ["template.0","template.1","template.2","template.3"];
my $templatedata = [
	["template.0","gram_positive_template","genome_scale_model","Bacteria","0","1","chenry"],
	["template.1","gram_negative_template","genome_scale_model","Bacteria","0","1","chenry"],
	["template.2","core_template","core_model","Bacteria","0","1","chenry"],
	["template.3","plant_template","genome_scale_model","Plant","0","1","seaver"],
];
#Printing template biomasses reactions
open($fh, ">", $directory."TemplateBiomasses.tsv");
my $columns = [qw(
	id
	name
	type
	other
	dna
	rna
	protein
	lipid
	cellwall
	cofactor
	energy
	template_id
	template_name
	template_modeltype
	template_domain
	template_version
	template_is_current
	template_owner
	compartment_ids
	compound_ids
	compound_data
)];
print $fh join("\t",@{$columns})."\n";
for (my $i=0; $i < @{$templates}; $i++) {
	my $biomasses = $templates->[$i]->templateBiomasses();
	for (my $j=0; $j < @{$biomasses}; $j++) {
		my $compounds = {};
		my $comps = {};
		my $bio = $biomasses->[$j];
		my $biocpds = $bio->templateBiomassComponents();
		my $compounddata = "";
		for (my $k=0; $k < @{$biocpds}; $k++) {
			my $biocpd = $biocpds->[$k];
			my $links = "";
			my $linkrefs = $biocpd->linked_compounds();
			for (my $m=0; $m < @{$linkrefs}; $m++) {
				if (length($links) > 0) {
					$links .= "|";
				}
				$links .= $linkrefs->[$m]->id()."{".$biocpd->link_coefficients()->[$m]."}";
			}
			if (length($compounddata) > 0) {
				$compounddata .= ";";
			}
			$compounds->{$biocpd->compound()->id()} = 1;
			$comps->{$biocpd->compartment()->id()} = 1;
			$compounddata .= $biocpd->compound()->id().":\"".$biocpd->compound()->name()."\":".$biocpd->coefficient().":".$biocpd->coefficientType().":".$biocpd->class().":".$links;
		}
		my $data = [
	    	$templatedata->[$i]->[0].".".$bio->id(),
	    	$bio->name(),
	    	"growth",
	    	$bio->other(),
	    	$bio->dna(),
	    	$bio->rna(),
	    	$bio->protein(),
	    	$bio->lipid(),
	    	$bio->cellwall(),
	    	$bio->cofactor(),
	    	$bio->energy(),
	    	$templatedata->[$i]->[0],
	    	$templatedata->[$i]->[1],
	    	$templatedata->[$i]->[2],
	    	$templatedata->[$i]->[3],
	    	$templatedata->[$i]->[4],
	    	$templatedata->[$i]->[5],
	    	$templatedata->[$i]->[6],
	    	"0:".join(";0:",keys(%{$comps})),
	    	join(";",keys(%{$compounds})),
	    	$compounddata
	    ];
	    print $fh join("\t",@{$data})."\n";
	}
}
close($fh);
#Printing complex roles
open($fh, ">", $directory."ComplexRoles.tsv");
$columns = [qw(
	complex_id
	complex_name
	complex_source
	complex_type
	role_id
	role_name
	role_type
	role_source
	role_aliases
	role_exemplar
	type
	triggering
	optional
)];
print $fh join("\t",@{$columns})."\n";
my $cpxs = $map->complexes();
my $idhash;
for (my $i=0; $i < @{$cpxs}; $i++) {
	my $cpx = $cpxs->[$i];
	$idhash->{$cpx->id()} = "cpx.".$i;
	my $cpxroles = $cpx->complexroles();
	for (my $j=0; $j < @{$cpxroles}; $j++) {
		my $cpxrole = $cpxroles->[$j];
		my $roleid = $cpxrole->role()->id();
		$roleid =~ s/ms//;
		my $data = [
	    	"cpx.".$i,
	    	"cpx.".$i,
	    	"ModelSEED",
	    	"SEED_role_complex",
	    	$roleid,
	    	$cpxrole->role()->name(),
	    	"SEED_role",
	    	"SEED",
	    	"searchname:".$cpxrole->role()->searchname(),
	    	"null",
	    	"role_mapping",
	    	$cpxrole->triggering(),
	    	$cpxrole->optionalRole(),
	    ];
	    print $fh join("\t",@{$data})."\n";
	}
}
close($fh);
#Printing compounds
my $cpdhash = {};
open($fh, ">", $directory."Compounds.tsv");
$columns = [qw(
	id
	abbreviation
	name
	formula
	mass
	source
	structure
	charge
	is_core
	is_obsolete
	linked_compound
	is_cofactor
	deltag
	deltagerr
	pka
	pkb
	abstract_compound
	comprised_of
	aliases
)];
print $fh join("\t",@{$columns})."\n";
my $cpds = $bio->compounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	my $cpd = $cpds->[$i];
	my $aliases = "";
	my $pka = "";
	my $pkb = "";
	my $structure = "null";
	if (defined($cpd_structure->{$cpd->id()})) {
		$structure = $cpd_structure->{$cpd->id()};
	}
	my $abstractcpd = "null";
	if (defined($cpd->abstractCompound_ref())) {
		$abstractcpd = $cpd->abstractCompound()->id();
	}
	foreach my $atom (keys(%{$cpd->pkas()})) {
		foreach my $value (@{$cpd->pkas()->{$atom}}) {
			if (length($pka) > 0) {
				$pka .= ";";
			}
			$pka .= $value.":".$atom;
		}
	}
	foreach my $atom (keys(%{$cpd->pkbs()})) {
		foreach my $value (@{$cpd->pkbs()->{$atom}}) {
			if (length($pkb) > 0) {
				$pkb .= ";";
			}
			$pkb .= $value.":".$atom;
		}
	}
	my $aliasehash = $cpd->parent()->compound_aliases()->{$cpd->id()};
    foreach my $type (keys(%{$aliasehash})) {
    	foreach my $alias (@{$aliasehash->{$type}}) {
    		if (length($aliases) > 0) {
	    		$aliases .= ";";
	    	}
    		$aliases .= "\"".$type.":".$alias."\"";;
    	}
    }
    my $data = [
    	$cpd->id(),
    	$cpd->abbreviation(),
    	$cpd->name(),
    	$cpd->formula(),
    	$cpd->mass(),
    	"ModelSEED",
    	$structure,
    	$cpd->defaultCharge(),
    	1,
    	0,
    	"null",
    	$cpd->isCofactor(),
    	$cpd->deltaG(),
    	$cpd->deltaGErr(),
    	$pka,
    	$pkb,
    	$abstractcpd,
    	"null",
    	$aliases
    ];
	print $fh join("\t",@{$data})."\n";
	$cpdhash->{$cpd->id()} = 1;
}
$cpds = $pbio->compounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	my $cpd = $cpds->[$i];
	if (!defined($cpdhash->{$cpd->id()})) {
		my $aliases = "";
		my $pka = "";
		my $pkb = "";
		my $structure = "null";
		if (defined($cpd_structure->{$cpd->id()})) {
			$structure = $cpd_structure->{$cpd->id()};
		}
		my $abstractcpd = "null";
		if (defined($cpd->abstractCompound_ref())) {
			$abstractcpd = $cpd->abstractCompound()->id();
		}
		foreach my $atom (keys(%{$cpd->pkas()})) {
			foreach my $value (@{$cpd->pkas()->{$atom}}) {
				if (length($pka) > 0) {
					$pka .= ";";
				}
				$pka .= $atom.":".$value;
			}
		}
		foreach my $atom (keys(%{$cpd->pkbs()})) {
			foreach my $value (@{$cpd->pkbs()->{$atom}}) {
				if (length($pkb) > 0) {
					$pkb .= ";";
				}
				$pkb .= $atom.":".$value;
			}
		}
		my $aliasehash = $cpd->parent()->compound_aliases()->{$cpd->id()};
	    foreach my $type (keys(%{$aliasehash})) {
	    	foreach my $alias (@{$aliasehash->{$type}}) {
	    		if (length($aliases) > 0) {
		    		$aliases .= ";";
		    	}
	    		$aliases .= "\"".$type.":".$alias."\"";;
	    	}
	    }
		my $data = [
	    	$cpd->id(),
	    	$cpd->abbreviation(),
	    	$cpd->name(),
	    	$cpd->formula(),
	    	$cpd->mass(),
	    	"ModelSEED",
	    	$structure,
	    	$cpd->defaultCharge(),
	    	1,
	    	0,
	    	"null",
	    	$cpd->isCofactor(),
	    	$cpd->deltaG(),
	    	$cpd->deltaGErr(),
	    	$pka,
	    	$pkb,
	    	$abstractcpd,
	    	"null",
	    	$aliases
	    ];
		print $fh join("\t",@{$data})."\n";
		$cpdhash->{$cpd->id()} = 1;
	}
}
close($fh);
#Printing reactions
my $rxnhash = {};
open($fh, ">", $directory."Reactions.tsv");
$columns = [qw(
	id
	abbreviation
	name
	code
	stoichiometry
	is_transport
	equation
	definition
	reversibility
	direction
	abstract_reaction
	pathways
	aliases
	ec_numbers
	deltag
	deltagerr
	compound_ids
)];
print $fh join("\t",@{$columns})."\n";
my $rxns = $bio->reactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	my $rxn = $rxns->[$i];
	my $aliases = "";
	my $pathways = "";
	if (defined($rxn_pathways->{$rxn->id()})) {
		foreach my $type (keys(%{$rxn_pathways->{$rxn->id()}})) {
			foreach my $path (keys(%{$rxn_pathways->{$rxn->id()}->{$type}})) {
				if (length($pathways) > 0) {
					$pathways .= ";";
				}
				$pathways .= $type.":".$path;
			}
		}
	}
	my $ecnums = "null";
	my $stoichiometry = "";
	my $compounds = {};
	my $rgts = $rxn->reagents();
	for (my $j=0; $j < @{$rgts}; $j++) {
		if (length($stoichiometry) > 0) {
			$stoichiometry .= ";";
		}
		$stoichiometry .= $rgts->[$j]->coefficient().":".$rgts->[$j]->compound()->id().":".$rgts->[$j]->compartment()->id().":0:\"".$rgts->[$j]->compound()->name()."\"";
		$compounds->{$rgts->[$j]->compound()->id()} = 1;
	}
	my $abstractrxn = "null";
	if (defined($rxn->abstractReaction_ref())) {
		$abstractrxn = $rxn->abstractReaction()->id();
	}
	my $aliasehash = $rxn->parent()->reaction_aliases()->{$rxn->id()};
    foreach my $type (keys(%{$aliasehash})) {
    	if ($type eq "Enzyme Class") {
    		for (my $m=0; $m < @{$aliasehash->{$type}}; $m++) {
    			if ($aliasehash->{$type}->[$m] =~ m/\d+\.\d+\.\d+\.\d+/) {
    				if ($ecnums eq "null") {
    					$ecnums = $aliasehash->{$type}->[$m];
    				} else {
    					$ecnums .= ";".$aliasehash->{$type}->[$m];
    				}
    			}
    		}
	    } else {
    		foreach my $alias (@{$aliasehash->{$type}}) {
    			if (length($aliases) > 0) {
	    			$aliases .= ";";
	    		}
    			$aliases .= "\"".$type.":".$alias."\"";;
    		}
    	}
    }
    my $data = [
    	$rxn->id(),
    	$rxn->abbreviation(),
    	$rxn->name(),
    	$rxn->code(),
    	$stoichiometry,
    	$rxn->isTransport(),
    	$rxn->equation(),
    	$rxn->definition(),
    	$rxn->thermoReversibility(),
    	$rxn->direction(),
    	$abstractrxn,
    	$pathways,
    	$aliases,
    	$ecnums,
    	$rxn->deltaG(),
    	$rxn->deltaGErr(),
    	join(";",keys(%{$compounds}))
    ];
    print $fh join("\t",@{$data})."\n";
    $rxnhash->{$rxn->id()} = 1;
}
close($fh);
#Printing template reactions
open($fh, ">", $directory."TemplateReactions.tsv");
$columns = [qw(
	id
	reaction_id
	abbreviation
	name
	code
	stoichiometry
	is_transport
	equation
	definition
	model_direction
	gapfill_direction
	type
	base_cost
	forward_penalty
	reverse_penalty
	pathways
	aliases
	ec_numbers
	deltag
	deltagerr
	template_id
	template_name
	template_modeltype
	template_domain
	template_version
	template_is_current
	template_owner
	compartment_ids
	complex_ids
	compound_ids
)];
print $fh join("\t",@{$columns})."\n";
for (my $i=0; $i < @{$templates}; $i++) {
	my $rxns = $templates->[$i]->templateReactions();
	for (my $j=0; $j < @{$rxns}; $j++) {
		my $rxn = $rxns->[$j];
		my $aliases = "";
		my $pathways = "";
		if (defined($rxn_pathways->{$rxn->reaction()->id()})) {
			foreach my $type (keys(%{$rxn_pathways->{$rxn->reaction()->id()}})) {
				foreach my $path (keys(%{$rxn_pathways->{$rxn->reaction()->id()}->{$type}})) {
					if (length($pathways) > 0) {
						$pathways .= ";";
					}
					$pathways .= $type.":".$path;
				}
			}
		}
		my $complexes = {};
		my $cpxs = $rxn->complexs();
		for (my $j=0; $j < @{$cpxs}; $j++) {
			$complexes->{$idhash->{$cpxs->[$j]->id()}} = 1;
		}
		my $compounds = {};
		my $comps = {};
		my $stoichiometry = "";
		my $rgts = $rxn->reaction()->reagents();
		for (my $j=0; $j < @{$rgts}; $j++) {
			if (length($stoichiometry) > 0) {
				$stoichiometry .= ";";
			}
			$stoichiometry .= $rgts->[$j]->coefficient().":".$rgts->[$j]->compound()->id().":".$rgts->[$j]->compartment()->id().":0:\"".$rgts->[$j]->compound()->name()."\"";
			$compounds->{$rgts->[$j]->compound()->id()} = 1;
			$comps->{$rgts->[$j]->compartment()->id()} = 1;
		}
		my $abstractrxn = "null";
		if (defined($rxn->reaction()->abstractReaction_ref())) {
			$abstractrxn = $rxn->reaction()->abstractReaction()->id();
		}
		my $aliasehash = $rxn->reaction()->parent()->reaction_aliases()->{$rxn->reaction()->id()};
	    my $ecnums = "null";
	    foreach my $type (keys(%{$aliasehash})) {
	    	if ($type eq "Enzyme Class") {
	    		for (my $m=0; $m < @{$aliasehash->{$type}}; $m++) {
	    			if ($aliasehash->{$type}->[$m] =~ m/\d+\.\d+\.\d+\.\d+/) {
	    				if ($ecnums eq "null") {
	    					$ecnums = $aliasehash->{$type}->[$m];
	    				} else {
	    					$ecnums .= ";".$aliasehash->{$type}->[$m];
	    				}
	    			}
	    		}
	    	} else {
	    		foreach my $alias (@{$aliasehash->{$type}}) {
	    			if (length($aliases) > 0) {
		    			$aliases .= ";";
		    		}
	    			$aliases .= "\"".$type.":".$alias."\"";
	    		}
	    	}
	    }
	    my $compid = "c0";
	    my $data = [
	    	$templatedata->[$i]->[0].".".$rxn->reaction()->id()."_".$compid,
	    	$rxn->reaction()->id(),
	    	$rxn->reaction()->abbreviation()."_".$compid,
	    	$rxn->reaction()->name()."_".$compid,
	    	$rxn->reaction()->code(),
	    	$stoichiometry,
	    	$rxn->reaction()->isTransport(),
	    	$rxn->reaction()->equation(),
	    	$rxn->reaction()->definition(),
	    	$rxn->direction(),
	    	$rxn->GapfillDirection(),
	    	"null",
	    	$rxn->base_cost(),
	    	$rxn->forward_penalty(),
	    	$rxn->reverse_penalty(),
	    	$pathways,
	    	$aliases,
	    	$ecnums,
	    	$rxn->reaction()->deltaG(),
	    	$rxn->reaction()->deltaGErr(),
	    	$templatedata->[$i]->[0],
	    	$templatedata->[$i]->[1],
	    	$templatedata->[$i]->[2],
	    	$templatedata->[$i]->[3],
	    	$templatedata->[$i]->[4],
	    	$templatedata->[$i]->[5],
	    	$templatedata->[$i]->[6],
	    	"0:".join(";0:",keys(%{$comps})),
	    	join(";",keys(%{$complexes})),
	    	join(";",keys(%{$compounds}))
	    ];
	    print $fh join("\t",@{$data})."\n";
	}
}
close($fh);