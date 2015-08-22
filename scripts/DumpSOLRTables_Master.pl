#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long::Descriptive;

my ($opt, $usage) = describe_options("%c %o <directory>",
	[ "compounds=s", "path to master compounds file", { default => "../Biochemistry/compounds.master.tsv" } ],
	[ "compartments=s", "path to master compartments file", { default => "../Biochemistry/compartments.master.tsv" } ],
	[ "reactions=s", "path to master reactions file", { default => "../Biochemistry/reactions.master.tsv" } ],
	[ "aliasdir=s", "path to directory with alias files", { default => "../Aliases/" } ],
	[ "pathwaydir=s", "path to directory with pathway files", { default => "../Pathways/" } ],
	[ "structuredir=s", "path to directory with structure files", { default => "../Structures/" } ],
	[ "master=s", "path to output master biochemistry json file", { default => "../Biochemistry/biochemistry.master.json" } ],
	[ "help|h", "print usage message and exit" ]
);
print($usage->text), exit if $opt->help;
my $directory = $ARGV[0];
exit if !$directory || !-d $directory;
$directory.="/" if $directory !~ /\/$/;

my @temp=();
my $header = 1;

#######################################################
#Initialization
#######################################################

#Collect Aliases
opendir(my $AliasDir, $opt->aliasdir);
my @Files = grep { $_ =~ /\.aliases$/ } readdir($AliasDir);
closedir($AliasDir);

my $rxn_alias_hash = {};
my $p_rxn_alias_hash = {};
my $cpd_alias_hash = {};
my $p_cpd_alias_hash = {};

my %alias_cpd_hash = ();
my %alias_rxn_hash = ();
foreach my $file (sort @Files){
    $file =~ /^(\w+)\.aliases/;
    my $aliasSet = $1;
    $aliasSet = join(" ", split(/_/,$aliasSet)) if $aliasSet eq "Enzyme_Class";

    open(FH, "< ".$opt->aliasdir.$file);
    $header = 1;
    while(<FH>){
	chomp;
	if($header){$header--;next}
	@temp=split(/\t/,$_,-1);

	if($temp[1] =~ /^cpd/ || $temp[2] =~ /^cpd/){
	    foreach my $cpd (split(/\|/,$temp[1])){
		$cpd_alias_hash->{$aliasSet}->{$temp[0]}->{$cpd}=1;
		$alias_cpd_hash{$cpd}{$aliasSet}{$temp[0]}=1;
	    }
	    foreach my $cpd (split(/\|/,$temp[2])){
		$p_cpd_alias_hash->{$aliasSet}->{$temp[0]}->{$cpd}=1;

		#Need to revise the decision to forgo aliases found in more recent database
		if(!exists($alias_cpd_hash{$cpd}) && !exists($alias_cpd_hash{$cpd}{$aliasSet})){
		    $alias_cpd_hash{$cpd}{$aliasSet}{$temp[0]}=1;
		}
	    }
	}
	if($temp[1] =~ /^rxn/ || $temp[2] =~ /^rxn/){
	    foreach my $rxn (split(/\|/,$temp[1])){
		$rxn_alias_hash->{$aliasSet}->{$temp[0]}->{$rxn}=1;
		$alias_rxn_hash{$rxn}{$aliasSet}{$temp[0]}=1;
	    }
	    foreach my $rxn (split(/\|/,$temp[2])){
		$p_rxn_alias_hash->{$aliasSet}->{$temp[0]}->{$rxn}=1;

		#Need to revise the decision to forgo aliases found in more recent database
		if(!exists($alias_rxn_hash{$rxn}) && !exists($alias_rxn_hash{$rxn}{$aliasSet})){
		    $alias_rxn_hash{$rxn}{$aliasSet}{$temp[0]}=1;
		}
	    }
	}

    }
    close(FH);
}

my $rxn_pathways = {};
open(my $fh, "< ".$opt->pathwaydir."HopeScenarios.txt");
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
open($fh, "< ".$opt->pathwaydir."plantdefault.pathways.tsv");
my @headers = split(/\t/,<$fh>);
shift(@headers);
chomp($headers[$#headers]);
while(<$fh>){
    chomp;
    @temp=split(/\t/,$_,-1);
    my $id = shift (@temp);
    for(my $i=0;$i<scalar(@headers);$i++){
	next if $temp[$i] eq "null";
	foreach my $path (split(/\|/,$temp[$i])){
	    $rxn_pathways->{$id}{$headers[$i]}{$temp[$i]}=1;
	}
    }
}
close($fh);

my $cpd_structure = {};
open($fh, "< ".$opt->structuredir."KEGG_Charged_InChI.txt");
while (my $line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	if (defined($cpd_alias_hash->{KEGG}->{$array->[0]})) {
	    foreach my $cpdid (keys(%{$cpd_alias_hash->{KEGG}->{$array->[0]}})){
		if (!defined($cpd_structure->{$cpdid})) {
		    $cpd_structure->{$cpdid} = $array->[1];
		}
	    } 
	}

	if (defined($p_cpd_alias_hash->{KEGG}->{$array->[0]})) {
	    foreach my $cpdid (keys(%{$p_cpd_alias_hash->{KEGG}->{$array->[0]}})){
		if (!defined($cpd_structure->{$cpdid})) {
		    $cpd_structure->{$cpdid} = $array->[1];
		}
	    }
	}
}
close($fh);

open($fh, "< ".$opt->structuredir."MetaCyc_Charged_InChI.txt");
while (my $line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	if (defined($cpd_alias_hash->{MetaCyc}->{$array->[0]})) {
	    foreach my $cpdid (keys(%{$cpd_alias_hash->{MetaCyc}->{$array->[0]}})){
		if (!defined($cpd_structure->{$cpdid})) {
		    $cpd_structure->{$cpdid} = $array->[1];
		}
	    } 
	}

	if (defined($p_cpd_alias_hash->{MetaCyc}->{$array->[0]})) {
	    foreach my $cpdid (keys(%{$p_cpd_alias_hash->{MetaCyc}->{$array->[0]}})){
		if (!defined($cpd_structure->{$cpdid})) {
		    $cpd_structure->{$cpdid} = $array->[1];
		}
	    }
	}
}
close($fh);

#Retreiving templates
# Need to get these from source file
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
#Printing complex roles
open($fh, ">", $directory."ComplexRoles.tsv");
my $columns = [qw(
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
my $cpxs = $map->complexes(); # Need to get these from source file
my $idhash;
my $count=0;
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
	$count = $i;
}

$cpxs = $pmap->complexes(); # Need to get these from source file
for (my $i=0; $i < @{$cpxs}; $i++) {
    $count++;
	my $cpx = $cpxs->[$i];
	$idhash->{$cpx->id()} = "cpx.".$count;
	my $cpxroles = $cpx->complexroles();
	for (my $j=0; $j < @{$cpxroles}; $j++) {
		my $cpxrole = $cpxroles->[$j];
		my $roleid = $cpxrole->role()->id();
		$roleid =~ s/ms//;
		my $data = [
	    	"cpx.".$count,
	    	"cpx.".$count,
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
#As it stands, it's a copy of the master compounds file with the aliases integrated
open(FH, "< ".$opt->compounds);
open($fh, ">", $directory."Compounds.tsv");
$header = 1;
undef(@headers);
my %Compounds=();
while(<FH>){
    chomp;
    if($header){
	@headers = split(/\t/,$_,-1);
	print $fh $_."\n";
	$header--;
	next;
    }
    @temp=split(/\t/,$_,-1);

    #map values to keys
    #probably not that necessary, but useful if column order changes
    my %cpdHash=();
    for(my $i=0;$i<scalar(@headers);$i++){
	$cpdHash{$headers[$i]}=$temp[$i];
    }

    my @aliases = ();
    foreach my $aliasSet (keys %{$alias_cpd_hash{$cpdHash{id}}}){
	foreach my $alias (keys %{$alias_cpd_hash{$cpdHash{id}}{$aliasSet}}){
	    push(@aliases, "\"".$aliasSet.":".$alias."\"");
	}
    }

    $cpdHash{aliases}= scalar(@aliases)>0 ? join(";",@aliases) : "null";

    print $fh join("\t", map { $cpdHash{$_} } @headers),"\n";

    $Compounds{$cpdHash{id}}=\%cpdHash;
}
close($fh);

#Printing reactions
#As it stands, it's a copy of the master reactions file with the pathways, aliases, and ec numbers integrated
open(FH, "< ".$opt->reactions);
open($fh, ">", $directory."Reactions.tsv");
$header = 1;
undef(@headers);
my %Reactions=();
while(<FH>){
    chomp;
    if($header){
	@headers = split(/\t/,$_,-1);
	print $fh  join("\t", grep { $_ ne 'is_obsolete' && $_ ne 'linked_reaction' } @headers),"\n";
	$header--;
	next;
    }
    @temp=split(/\t/,$_,-1);

    #map values to keys
    #probably not that necessary, but useful if column order changes
    my %rxnHash=();
    for(my $i=0;$i<scalar(@headers);$i++){
	$rxnHash{$headers[$i]}=$temp[$i];
    }

    my @ecnums = ();		
    my @aliases = ();
    foreach my $aliasSet (keys %{$alias_rxn_hash{$rxnHash{id}}}){
	foreach my $alias (keys %{$alias_rxn_hash{$rxnHash{id}}{$aliasSet}}){
	    #Only include full ec numbers (?)
	    if ($aliasSet eq "Enzyme Class"){
		if($alias =~ m/\d+\.\d+\.\d+\.\d+/){
		    push(@ecnums, $alias);
		}
	    }else{
		push(@aliases, "\"".$aliasSet.":".$alias."\"");
	    }
	}
    }

    $rxnHash{aliases}= scalar(@aliases)>0 ? join(";",@aliases) : "null";
    $rxnHash{ec_numbers}= scalar(@ecnums)>0 ? join(";",@ecnums) : "null";

    my @pathways = ();
    if (defined($rxn_pathways->{$rxnHash{id}})) {
	foreach my $type (keys(%{$rxn_pathways->{$rxnHash{id}}})) {
	    foreach my $path (keys(%{$rxn_pathways->{$rxnHash{id}}{$type}})) {
		push(@pathways, $type.":".$path);
	    }
	}
    }

    $rxnHash{pathways}= scalar(@pathways)>0 ? join(";",@pathways) : "null";

    print $fh join("\t", map { $rxnHash{$_} } grep { $_ ne 'is_obsolete' && $_ ne 'linked_reaction' } @headers),"\n";
    $Reactions{$rxnHash{id}}=\%rxnHash;
}
close($fh);
#Printing template biomasses reactions
open($fh, ">", $directory."TemplateBiomasses.tsv");
$columns = [qw(
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
		my @compounddata = ();
		for (my $k=0; $k < @{$biocpds}; $k++) {
			my $biocpd = $biocpds->[$k];
			my $biocpd_id = $biocpd->compound()->id();
			my @links = ();
			my $linkrefs = $biocpd->linked_compounds();
			for (my $m=0; $m < @{$linkrefs}; $m++) {
			    push(@links, $linkrefs->[$m]->id()."{".$biocpd->link_coefficients()->[$m]."}");
			}
			$compounds->{$biocpd_id} = 1;
			$comps->{$biocpd->compartment()->id()} = 1;
			push(@compounddata, $biocpd_id.":\"".$Compounds{$biocpd_id}{name}."\":".$biocpd->coefficient().":".$biocpd->coefficientType().":".$biocpd->class().":".join("|",@links));
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
		join(";",@compounddata)
	    ];
	    print $fh join("\t",@{$data})."\n";
	}
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
		my $complexes = {};
		my $cpxs = $rxn->complexs();
		for (my $j=0; $j < @{$cpxs}; $j++) {
			$complexes->{$idhash->{$cpxs->[$j]->id()}} = 1;
		}
		my $compounds = {};
		my $comps = {};
#		my $rgts = [split(/;/,$Reactions{$rxn->reaction()->id()}{stoichiometry})];
		my $rgts = $rxn->reaction()->reagents();
		for (my $j=0; $j < @{$rgts}; $j++) {
#		    my ($coef,$cpd,$cmpt) = split(/:/,$rgts->[$j]);
		    my ($cpd,$cmpt) = ($rgts->[$j]->compound()->id(),$rgts->[$j]->compartment()->id());
		    $compounds->{$cpd}=1;
		    $comps->{$cmpt}=1;
		}

		my $rxn_id = $rxn->reaction()->id();
	    my $compid = "c0";
	    my $data = [
	    	$templatedata->[$i]->[0].".".$rxn_id."_".$compid,
	    	$rxn_id,
	    	$Reactions{$rxn_id}{abbreviation}."_".$compid,
	    	$Reactions{$rxn_id}{name}."_".$compid,
#	    	$Reactions{$rxn_id}{code},
#	    	$Reactions{$rxn_id}{stoichiometry},
		$rxn->reaction()->code(),
		$rxn->reaction()->stoichiometry(),
	    	$Reactions{$rxn_id}{is_transport},
#	    	$Reactions{$rxn_id}{equation},
#	    	$Reactions{$rxn_id}{definition},
		$rxn->reaction()->equation(),
		$rxn->reaction()->definition(),
	    	$rxn->direction(),
	    	defined($rxn->GapfillDirection()) ? $rxn->GapfillDirection() : "=",
	    	"null",
	    	defined($rxn->base_cost()) ? $rxn->base_cost() : 0,
	    	defined($rxn->forward_penalty()) ? $rxn->base_cost() : 0,
	    	defined($rxn->reverse_penalty()) ? $rxn->base_cost() : 0,
	    	$Reactions{$rxn_id}{pathways},
	    	$Reactions{$rxn_id}{aliases},
	    	$Reactions{$rxn_id}{ec_numbers},
		$Reactions{$rxn_id}{deltag},
		$Reactions{$rxn_id}{deltagerr},
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
