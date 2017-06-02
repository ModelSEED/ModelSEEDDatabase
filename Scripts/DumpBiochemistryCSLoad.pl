#!/usr/bin/perl -w

use strict;
use Config::Simple;
use Bio::KBase::workspaceService::Client;
use Bio::KBase::fbaModelServices::Impl;

my $config = $ARGV[0];
my $directory = $ARGV[1];
if (!defined($config)) {
	print STDERR "No config file provided!\n";
	exit(-1);
}
if (!-e $config) {
	print STDERR "Config file ".$config." not found!\n";
	exit(-1);
}
#Params: writesbml.wsurl, writesbml.fbaurl, writesbml.auth
my $c = Config::Simple->new();
$c->read($config);
my $wserv = Bio::KBase::workspaceService::Client->new($c->param("kbclientconfig.wsurl"));
my $fba = Bio::KBase::fbaModelServices::Impl->new({
	"fba-url" => "",
	"probanno-url" => "",
	"mssserver-url" => "http://bio-data-1.mcs.anl.gov/services/ms_fba",
	accounttype => "kbase",
	"workspace-url" => $c->param("kbclientconfig.wsurl"),
	defaultJobState => "queued",
	"idserver-url" => "http://bio-data-1.mcs.anl.gov/services/idserver"
});
$fba->_setContext(undef,{auth => $c->param("kbclientconfig.auth")});
my $bio = $fba->_get_msobject("Biochemistry","PlantSEED","KBase-Biochem");
my $map = $fba->_get_msobject("Mapping","kbase","default");
my $pmap = $fba->_get_msobject("Mapping","PlantSEED","PlantSEED");
#    LOCATION ( a.k.a. COMPARTMENT )
#    Location: A location is a place where a reaction's compounds can originate or
#              end up (e.g. cell wall, extracellular, cytoplasm).
#
#        Table: Location
#            id (string): Unique identifier for this Location.
#            hierarchy (int): a number indicating where this location occurs in relation
#                             other locations in the cell. Zero indicates extra-cellular.
#            mod-date (date): date and time of the last modification to the compartment's
#                             definition
#            msid (string): common modeling ID of this location
#            name (string): common name for the location
open(my $fh, ">", $directory."location.dtx");
print $fh "id\thierarchy\tmod-date\tmsid\tname\n";
my $locationID = {
	c => "kb|loc.3",
	v => "kb|loc.11",
	p => "kb|loc.2",
	n => "kb|loc.7",
	h => "kb|loc.8",
	d => "kb|loc.12",
	l => "kb|loc.6",
	e => "kb|loc.0",
	m => "kb|loc.9",
	r => "kb|loc.5",
	x => "kb|loc.10",
	w => "kb|loc.1",
	g => "kb|loc.4"
};
my $comps = $bio->compartments();
for (my $i=0; $i < @{$comps}; $i++) {
	my $comp = $comps->[$i];
	print $fh $locationID->{$comp->id()}."\t".$comp->hierarchy()."\t".$comp->modDate()."\t".$comp->id()."\t".$comp->name()."\n";
}
close($fh);
#    Complex: A complex is a set of chemical reactions that act in concert to effect
#             a role.
#    
#        Table: Complex
#            id (string): Unique identifier for this Complex.
#            mod-date (date): date and time of the last change to this complex's
#                             definition
my $cpxhash;
my $cpxidhash;
my $mergedcpx;
my $eliminatedcpx;
my $cpxname;
my $cpxs = $map->complexes();
for (my $i=0; $i < @{$cpxs}; $i++) {
	my $cpx = $cpxs->[$i];
	my $roles = $cpx->complexroles();
	my $rolelist = [];
	for (my $j=0; $j < @{$roles}; $j++) {
		push(@{$rolelist},$roles->[$j]->role()->name()."___".$roles->[$j]->optional()."___".$roles->[$j]->type()."___".$roles->[$j]->triggering());
	}
	my $uniqueid = join(";",sort(@{$rolelist}));
	if (length($uniqueid) == 0) {
		$eliminatedcpx ->{$cpx->name()} = 1;
	} else {
		if (defined($cpxhash->{$uniqueid})) {
			print "Merged:".$cpx->name()."|".$cpxhash->{$uniqueid}->name()."\n";
			$mergedcpx->{$cpx->name()} = $cpxhash->{$uniqueid}->name();
		} else {
			$cpxname->{$uniqueid} = $roles->[0]->role()->name();
			$cpxhash->{$uniqueid} = $cpx;
		}
		if (defined($cpxidhash->{$cpx->name()})) {
			print "ID collision:".$cpx->name()."\n";
		} else {
			$cpxidhash->{$cpx->name()} = $cpx;
		}
	}
}
$cpxs = $pmap->complexes();
for (my $i=0; $i < @{$cpxs}; $i++) {
	my $cpx = $cpxs->[$i];
	my $roles = $cpx->complexroles();
	my $rolelist = [];
	for (my $j=0; $j < @{$roles}; $j++) {
		push(@{$rolelist},$roles->[$j]->role()->name()."___".$roles->[$j]->optional()."___".$roles->[$j]->type()."___".$roles->[$j]->triggering());
	}
	my $uniqueid = join(";",sort(@{$rolelist}));
	if (length($uniqueid) == 0) {
		$eliminatedcpx ->{$cpx->name()} = 1;
	} else {
		if (defined($cpxhash->{$uniqueid})) {
			print "Merged:".$cpx->name()."|".$cpxhash->{$uniqueid}->name()."\n";
			$mergedcpx->{$cpx->name()} = $cpxhash->{$uniqueid}->name();
		} else {
			$cpxname->{$uniqueid} = $roles->[0]->role()->name();
			$cpxhash->{$uniqueid} = $cpx;
		}
		if (defined($cpxidhash->{$cpx->name()})) {
			print "ID collision:".$cpx->name()."\n";
		} else {
			$cpxidhash->{$cpx->name()} = $cpx;
		}
	}
}
open($fh, ">", $directory."complex.dtx");
print $fh "id\tmod-date\tsource-id\n";
foreach my $cpx (keys(%{$cpxhash})) {
	my $id = "kb|".$cpxhash->{$cpx}->name();
	$id =~ s/cpx0*/cpx./;
	$id =~ s/mscpx/cpx/;
	print $fh $id."\t".$cpxhash->{$cpx}->modDate()."\t".$cpxhash->{$cpx}->name()."\n";
}
close($fh);
#        Table: ComplexName
#            id (string): Unique identifier for this Complex.
#            name (string): name of this complex. Not all complexes have names.
open($fh, ">", $directory."complexName.dtx");
print $fh "id\tmsid\tname\n";
foreach my $cpx (keys(%{$cpxhash})) {
	my $id = "kb|".$cpxhash->{$cpx}->name();
	$id =~ s/cpx0*/cpx./;
	$id =~ s/mscpx/cpx/;
	print $fh $id."\t".$cpxname->{$cpx}."\n";
}
close($fh);
#    Compound: A compound is a chemical that participates in a reaction. Both ligands
#              and reaction components are treated as compounds.
#    
#        Table: Compound
#            id (string): Unique identifier for this Compound.
#            mass (float): atomic mass of the compound
#            mod-date (date): date and time of the last modification to the compound
#                             definition
#            ubiquitous (boolean): TRUE if this compound is found in most reactions,
#                                  else FALSE
#            abbr (string): shortened abbreviation for the compound name
#            default-charge (float): computed charge of the compound in a pH-neutral
#                                    solution
#            deltaG (float): Gibbs free-energy coefficient for this compound
#            deltaG-error (float): error bounds on the deltaG value
#            formula (string): a pH-neutral formula for the compound
#            label (string): primary name of the compound, for use in displaying
#                            reactions
#            msid (string): common modeling ID of this compound
#            uncharged-formula (string): a electrically neutral formula for the compound   
open($fh, ">", $directory."compound.dtx");
print $fh "abbr\tcharge\tdeltaG\tdeltaG-error\tformula\tid\tlabel\tmass\tmod-date\tsource-id\n";
my $cpds = $bio->compounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	my $cpd = $cpds->[$i];
	my $id = "kb|".$cpd->id();
	$id =~ s/cpd0*/cpd./;
	print $fh $cpd->abbreviation()."\t".$cpd->defaultCharge()."\t".$cpd->deltaG()."\t".$cpd->deltaGErr()."\t".$cpd->formula()."\t".$id."\t".$cpd->name()."\t".$cpd->mass()."\t".$cpd->modDate()."\t".$cpd->id()."\n";
}
close($fh);
#    Media: A media describes the chemical content of the solution in which cells
#           are grown in an experiment or for the purposes of a model. The key is
#           the common media name. The nature of the media is described by its relationship
#           to its constituent compounds.
#
#        Table: Media
#            id (string): Unique identifier for this Media.
#            is-minimal (boolean): TRUE if this media condition is considered minimal
#            mod-date (date): date and time of the last modification to the media's
#                             definition
#            name (string): descriptive name of the media
#            type (string): type of the medium (aerobic or anaerobic)
open($fh, ">", $directory."media.dtx");
print $fh "id\tis-minimal\tmod-date\tsource-id\tname\ttype\n";
my $medias = $bio->media();
for (my $i=0; $i < @{$medias}; $i++) {
	my $media = $medias->[$i];
	print $fh "kb|med.".($i+1)."\t".$media->isMinimal()."\t".$media->modDate()."\t".$media->name()."\t".$media->name()."\t".$media->type()."\n";
}
close($fh);
#    Reaction: A reaction is a chemical process that converts one set of compounds
#              (substrate) to another set (products).
#
#        Table: Reaction
#            id (string): Unique identifier for this Reaction.
#            default-protons (float): number of protons absorbed by this reaction
#                                     in a pH-neutral environment
#            deltaG (float): Gibbs free-energy coefficient for the reaction
#            deltaG-error (float): error bounds on the deltaG value
#            direction (char): direction of this reaction (> for forward-only, <
#                              for backward-only, = for bidirectional)
#            mod-date (date): date and time of the last modification to this reaction's
#                             definition
#            thermodynamic-reversibility (char): computed reversibility of this reaction
#                                                in a pH-neutral environment
#            abbr (string): abbreviated name of this reaction
#            msid (string): common modeling ID of this reaction
#            name (string): descriptive name of this reaction
#            status (string): string indicating additional information about this
#                             reaction, generally indicating whether the reaction
#                             is balanced and/or accurate
open($fh, ">", $directory."reaction.dtx");
print $fh "abbr\tdefault-protons\tdeltaG\tdeltaG-error\tdirection\tmod-date\tname\tid\tsource-id\tstatus\tthermodynamic-reversibility\n";
my $rxns = $bio->reactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	my $rxn = $rxns->[$i];
	my $id = "kb|".$rxn->id();
	$id =~ s/rxn0*/rxn./;
	print $fh $rxn->abbreviation()."\t".$rxn->defaultProtons()."\t".$rxn->deltaG()."\t".$rxn->deltaGErr()."\t".$rxn->direction()."\t".$rxn->modDate()."\t".$rxn->name()."\t".$id."\t".$rxn->id()."\t".$rxn->status()."\t".$rxn->thermoReversibility()."\n";
}
close($fh);
#    HasCompoundAliasFrom: This relationship connects a source (database or organization)
#                          with the compounds for which it has assigned names (aliases).
#                          The alias itself is stored as intersection data.
#    
#        Table: HasCompoundAliasFrom
#            from-link (string): id of the source Source.
#            to-link (string): id of the target Compound.
#            alias (string): alias for the compound assigned by the source
open($fh, ">", $directory."hasCompoundAliasFrom.dtx");
print $fh "alias\tto-link\tfrom-link\n";
my $cpdaliases = $bio->queryObjects("aliasSets", { attribute => "compounds" });
for (my $i=0; $i < @{$cpdaliases}; $i++) {
	my $aliases = $cpdaliases->[$i]->aliases;
	foreach my $als_name (keys %$aliases) {
		my $uuids = $aliases->{$als_name};
		foreach my $uuid (@$uuids) {
			my $id = "kb|".$bio->getObject("compounds",$uuid)->id();
			$id =~ s/cpd0*/cpd./;
			print $fh $als_name."\t".$id."\t".$cpdaliases->[$i]->name."\n";
		}
	}
}
close($fh);
#    HasPresenceOf: This relationship connects a media to the compounds that occur
#                   in it. The intersection data describes how much of each compound
#                   can be found.
#    
#        Table: HasPresenceOf
#            from-link (string): id of the source Media.
#            to-link (string): id of the target Compound.
#            concentration (float): concentration of the compound in the media
#            maximum-flux (float): maximum flux of the compound for this media
#            minimum-flux (float): minimum flux of the compound for this media    
open($fh, ">", $directory."hasPresenceOf.dtx");
print $fh "concentration\tfrom-link\tmaximum-flux\tminimum-flux\tto-link\n";
$medias = $bio->media();
for (my $i=0; $i < @{$medias}; $i++) {
	my $media = $medias->[$i];
	my $mediacpds = $media->mediacompounds();
	for (my $j=0; $j < @{$mediacpds}; $j++) {
		my $id = "kb|".$mediacpds->[$j]->compound()->id();
		$id =~ s/cpd0*/cpd./;
		print $fh $mediacpds->[$j]->concentration()."\t"."kb|med.".($i+1)."\t".$mediacpds->[$j]->maxFlux()."\t".$mediacpds->[$j]->minFlux()."\t".$id."\n";
	}
}
close($fh); 
#    HasReactionAliasFrom: This relationship connects a source (database or organization)
#                          with the reactions for which it has assigned names (aliases).
#                          The alias itself is stored as intersection data.
#    
#        Table: HasReactionAliasFrom
#            from-link (string): id of the source Source.
#            to-link (string): id of the target Reaction.
#            alias (string): alias for the reaction assigned by the source    
open($fh, ">", $directory."hasReactionAliasFrom.dtx");
print $fh "alias\tto-link\tfrom-link\n";
my $rxnaliases = $bio->queryObjects("aliasSets", { attribute => "reactions" });
for (my $i=0; $i < @{$rxnaliases}; $i++) {
	my $aliases = $rxnaliases->[$i]->aliases;
	foreach my $als_name (keys %$aliases) {
		my $uuids = $aliases->{$als_name};
		foreach my $uuid (@$uuids) {
			my $id = "kb|".$bio->getObject("reactions",$uuid)->id();
			$id =~ s/rxn0*/rxn./;
			print $fh $als_name."\t".$id."\t".$rxnaliases->[$i]->name."\n";
		}
	}
}
close($fh);
#    HasStep: This relationship connects a complex to the reaction instances that
#             work together to make the complex happen.
#    
#        Table: HasStep
#            from-link (string): id of the source Complex.
#            to-link (string): id of the target Reaction.
open($fh, ">", $directory."hasStep.dtx");
print $fh "from-link\tto-link\n";
my $temps = [
	$fba->_get_msobject("ModelTemplate","PlantSEED","PlantTemplate"),
	$fba->_get_msobject("ModelTemplate","KBaseTemplateModels","GramNegModelTemplate"),
	$fba->_get_msobject("ModelTemplate","KBaseTemplateModels","GramPosModelTemplate")
];
my $cpxrxnHash = {};
foreach my $temp (@{$temps}) {
	my $tmprxns = $temp->templateReactions();
	foreach my $tmprxn (@{$tmprxns}) {
		foreach my $cpx (@{$tmprxn->complexes()}) {
			if (!defined($eliminatedcpx->{$cpx->name()})) {
				if (defined($mergedcpx->{$cpx->name()})) {
					$cpxrxnHash->{$mergedcpx->{$cpx->name()}}->{$tmprxn->reaction()->id()} = 1;
				} else{
					$cpxrxnHash->{$cpx->name()}->{$tmprxn->reaction()->id()} = 1;
				}	
			}
		}
	}
}
foreach my $cpxrxn (keys(%{$cpxrxnHash})) {
	foreach my $rxn (keys(%{$cpxrxnHash->{$cpxrxn}})) {
		my $cpxrxn = "kb|".$cpxrxn;
		$cpxrxn =~ s/cpx0*/cpx./;
		$cpxrxn =~ s/mscpx/cpx/;
		my $id = "kb|".$rxn;
		$id =~ s/rxn0*/rxn./;
		print $fh $cpxrxn."\t".$id."\n";
	}
}
close($fh);
#    LocalizedCompound: This entity represents a compound occurring in a specific
#                       location. A reaction always involves localized compounds.
#                       If a reaction occurs entirely in a single location, it will
#                       frequently only be represented by the cytoplasmic versions
#                       of the compounds; however, a transport always uses specifically
#                       located compounds.
#
#        Table: LocalizedCompound
#            id (string): Unique identifier for this LocalizedCompound.
open($fh, ">", $directory."LocalizedCompound.dtx");
print $fh "id\n";
my $loccomphash = {};
$rxns = $bio->reactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	my $rxn = $rxns->[$i];
	my $cpds = $rxns->[$i]->reagents();
	for (my $j=0; $j < @{$cpds}; $j++) {
		my $id = $locationID->{$cpds->[$j]->compartment()->id()}.".".$cpds->[$j]->compound()->id();
		$id =~ s/cpd0*/cpd./;
		$loccomphash->{$id} = 1;
	}
}
foreach my $loccomp (keys(%{$loccomphash})) {
	print $fh $loccomp."\n";
}
close($fh);
#    INVOLVES ( a.k.a. REAGENT )
#    Involves: This relationship connects a reaction to the specific localized compounds
#              that participate in it.
#
#        Table: Involves
#            from-link (string): id of the source Reaction.
#            to-link (string): id of the target LocalizedCompound.
#            coefficient (float): Number of molecules of the compound that participate
#                                 in a single instance of the reaction. For example,
#                                 if a reaction produces two water molecules, the
#                                 stoichiometry of water for the reaction would be
#                                 two. When a reaction is written on paper in chemical
#                                 notation, the stoichiometry is the number next
#                                 to the chemical formula of the compound. The value
#                                 is negative for substrates and positive for products.
#            cofactor (boolean): TRUE if the compound is a cofactor; FALSE if it
#                                is a major component of the reaction.
#            is-transport (boolean): TRUE if the compound is being transported out
#                                    of or into the reaction's compartment, FALSE
#                                    if it stays in the same compartment
open($fh, ">", $directory."Involves.dtx");
print $fh "coefficient\tcofactor\tfrom-link\tto-link\n";
$rxns = $bio->reactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	my $rxn = $rxns->[$i];
	my $cpds = $rxns->[$i]->reagents();
	for (my $j=0; $j < @{$cpds}; $j++) {
		my $id = "kb|".$rxn->id();
		$id =~ s/rxn0*/rxn./;
		my $lid = $locationID->{$cpds->[$j]->compartment()->id()}.".".$cpds->[$j]->compound()->id();
		$lid =~ s/cpd0*/cpd./;
		print $fh $cpds->[$j]->coefficient()."\t".$cpds->[$j]->isCofactor()."\t".$id."\t".$lid."\n";
	}
}
close($fh);
#    IsTriggeredBy: A complex can be triggered by many roles. A role can trigger
#                   many complexes.
#    
#        Table: IsTriggeredBy
#            from-link (string): id of the source Complex.
#            to-link (string): id of the target Role.
#            optional (boolean): TRUE if the role is not necessarily required to
#                                trigger the complex, else FALSE
#            type (char): ask Chris
#            triggering (boolean): TRUE if the presence of the role requires including
#            the complex in the model, else FALSE.
open($fh, ">", $directory."isTriggeredBy.dtx");
print $fh "from-link\toptional\tto-link\ttriggering\ttype\n";
foreach my $cpx (keys(%{$cpxidhash})) {
	my $roles = $cpxidhash->{$cpx}->complexroles();
	foreach my $role (@{$roles}) {
		my $id = "kb|".$cpxidhash->{$cpx}->name();
		$id =~ s/cpx0*/cpx./;
		$id =~ s/mscpx/cpx/;
		print $fh $id."\t".$role->optional()."\t".$role->role()->name()."\t".$role->triggering()."\t".$role->type()."\n";
	}
}
close($fh);
#   IsParticipatingAt: This relationship connects a localized compound to the location
#                      in which it occurs during one or more reactions.
#
#       Table: IsParticipatingAt
#           from-link (string): id of the source Location.
#           to-link (string): id of the target LocalizedCompound.
#    
#    ParticipatesAs: This relationship connects a compound to the reagents that represent
#                    its participation in reactions.
#    
#        Table: ParticipatesAs
#            from-link (string): id of the source Compound.
#            to-link (string): id of the target LocalizedCompound.
open($fh, ">", $directory."IsParticipatingAt.dtx");
print $fh "from-link\tto-link\n";
$rxns = $bio->reactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	my $rxn = $rxns->[$i];
	my $cpds = $rxns->[$i]->reagents();
	for (my $j=0; $j < @{$cpds}; $j++) {
		my $id = $locationID->{$cpds->[$j]->compartment()->id()}.".".$cpds->[$j]->compound()->id();
		$id =~ s/cpd0*/cpd./;
		print $fh $locationID->{$cpds->[$j]->compartment()->id()}."\t".$id."\n";
	}
}
close($fh);
open($fh, ">", $directory."ParticipatesAs.dtx");
print $fh "from-link\tto-link\n";
$rxns = $bio->reactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	my $rxn = $rxns->[$i];
	my $cpds = $rxns->[$i]->reagents();
	for (my $j=0; $j < @{$cpds}; $j++) {
		my $id = $cpds->[$j]->compound()->id();
		$id =~ s/cpd0*/cpd./;
		my $lid = $locationID->{$cpds->[$j]->compartment()->id()}.".".$cpds->[$j]->compound()->id();
		$lid =~ s/cpd0*/cpd./;
		print $fh "kb|".$id."\t".$lid."\n";
	}
}
close($fh);
1;
