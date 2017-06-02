#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use JSON::XS;

my $arrayfields = {
	pka => 1,
	pkb => 1,
	linked_reaction => 1,
	compound_ids => 1,
	stoichiometry => 1,
};

my $keylist = [qw(
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
)];
my $compoundhash;
my $filename = $Bin."/../Biochemistry/compounds.master.tsv";
open(my $fh, "<", $filename);
my $heading = <$fh>;
chomp($heading);
my $array = [split(/\t/,$heading)];
my $headinghash;
for (my $i=0; $i <@{$array}; $i++) {
	$headinghash->{$array->[$i]} = $i;
}
while (my $line = <$fh>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	my $current = {};
	for (my $i=0; $i < @{$keylist}; $i++) {
		if (defined($headinghash->{$keylist->[$i]})) {
			if (defined($array->[$headinghash->{$keylist->[$i]}])) {
				if (length($array->[$headinghash->{$keylist->[$i]}]) > 0 && $array->[$headinghash->{$keylist->[$i]}] ne "null") {
					$current->{$keylist->[$i]} = $array->[$headinghash->{$keylist->[$i]}];
					if (defined($arrayfields->{$keylist->[$i]})) {
						$current->{$keylist->[$i]} = [split(/;/,$current->{$keylist->[$i]})];
					}
				}
			}
		}
	}
	$compoundhash->{$current->{id}} = $current;
}
close($fh);

$keylist = [qw(
id
abbreviation
name
code
stoichiometry
equation
definition
reversibility
direction
deltag
deltagerr
compound_ids
status
is_obsolete
linked_reaction
)];
my $kegghash;
my $metacychash;
my $reactionhash;
$filename = $Bin."/../Biochemistry/reactions.master.tsv";
open(my $fhh, "<", $filename);
$heading = <$fhh>;
chomp($heading);
$array = [split(/\t/,$heading)];
$headinghash = {};
for (my $i=0; $i <@{$array}; $i++) {
	$headinghash->{$array->[$i]} = $i;
}
while (my $line = <$fhh>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	my $current = {};
	for (my $i=0; $i < @{$keylist}; $i++) {
		if (defined($headinghash->{$keylist->[$i]})) {
			if (defined($array->[$headinghash->{$keylist->[$i]}])) {
				if (length($array->[$headinghash->{$keylist->[$i]}]) > 0 && $array->[$headinghash->{$keylist->[$i]}] ne "null") {
					$current->{$keylist->[$i]} = $array->[$headinghash->{$keylist->[$i]}];
					if (defined($arrayfields->{$keylist->[$i]})) {
						$current->{$keylist->[$i]} = [split(/;/,$current->{$keylist->[$i]})];
					}
				}
			}
		}
	}
	$reactionhash->{$current->{id}} = $current;
}
close($fh);

my $aliaslist = {
	kegg_aliases => "KEGG.aliases",
	names => "name.aliases",
	searchnames => "searchname.aliases",
	metacyc_aliases => "MetaCyc.aliases",
	bigg_aliases => "BiGG.aliases",
	ec_numbers => "Enzyme_Class.aliases"
};
foreach my $key (keys(%{$aliaslist})) {
	$filename = $Bin."/../Aliases/".$aliaslist->{$key};
	open(my $fhhh, "<", $filename);
	my $heading = <$fhhh>;
	chomp($heading);
	while (my $line = <$fhhh>) {
		chomp($line);
		$array = [split(/\t/,$line)];
		if (defined($array->[1]) && length($array->[1]) > 0) {
			if ($key eq "kegg_aliases") {
				if (!defined($kegghash->{$array->[0]})) {
					$kegghash->{$array->[0]} = $array->[1];
				}
			} elsif ($key eq "metacyc_aliases") {
				if (!defined($metacychash->{$array->[0]})) {
					$metacychash->{$array->[0]} = $array->[1];
				}
			}
			if (defined($compoundhash->{$array->[1]})) {
				push(@{$compoundhash->{$array->[1]}->{$key}},$array->[0]);
			} elsif (defined($reactionhash->{$array->[1]})) {
				push(@{$reactionhash->{$array->[1]}->{$key}},$array->[0]);
			}
		} elsif (defined($array->[2]) && length($array->[2]) > 0) {
			if ($key eq "kegg_aliases") {
				if (!defined($kegghash->{$array->[0]})) {
					$kegghash->{$array->[0]} = $array->[2];
				}
			} elsif ($key eq "metacyc_aliases") {
				if (!defined($metacychash->{$array->[0]})) {
					$metacychash->{$array->[0]} = $array->[2];
				}
			}
			if (defined($compoundhash->{$array->[2]})) {
				push(@{$compoundhash->{$array->[2]}->{$key}},$array->[0]);
			} elsif (defined($reactionhash->{$array->[2]})) {
				push(@{$reactionhash->{$array->[2]}->{$key}},$array->[0]);
			}
		}
	}
	close($fhhh);
}

$filename = $Bin."/../Structures/KEGG_Search_InChI.txt";
open(my $fa, "<", $filename);
while (my $line = <$fa>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	if (defined($array->[1]) && length($array->[1]) > 0) {
		if (defined($kegghash->{$array->[0]}) && defined($compoundhash->{$kegghash->{$array->[0]}})) {
			$compoundhash->{$kegghash->{$array->[0]}}->{search_inchi} = $array->[1];
		}
	}
}
close($fa);
$filename = $Bin."/../Structures/KEGG_Charged_InChI.txt";
open(my $fb, "<", $filename);
while (my $line = <$fb>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	if (defined($array->[1]) && length($array->[1]) > 0) {
		if (defined($kegghash->{$array->[0]}) && defined($compoundhash->{$kegghash->{$array->[0]}})) {
			$compoundhash->{$kegghash->{$array->[0]}}->{structure} = $array->[1];
		}
	}
}
close($fb);
$filename = $Bin."/../Structures/KEGG_Charged_MolAnalysis.tbl";
open(my $fc, "<", $filename);
while (my $line = <$fc>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	if (defined($array->[1]) && length($array->[1]) > 0) {
		if (defined($kegghash->{$array->[0]}) && defined($compoundhash->{$kegghash->{$array->[0]}})) {
			$compoundhash->{$kegghash->{$array->[0]}}->{groups} = [split(/\|/,$array->[2])];
			$compoundhash->{$kegghash->{$array->[0]}}->{charge} = $array->[3];
			$compoundhash->{$kegghash->{$array->[0]}}->{formula} = $array->[4];
			$compoundhash->{$kegghash->{$array->[0]}}->{mass} = $array->[6];
			$compoundhash->{$kegghash->{$array->[0]}}->{deltag} = $array->[7];
			$compoundhash->{$kegghash->{$array->[0]}}->{deltagerr} = $array->[8];
		}
	}
}
close($fc);

$filename = $Bin."/../Structures/MetaCyc_Search_InChI.txt";
open(my $fd, "<", $filename);
while (my $line = <$fd>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	if (defined($array->[1]) && length($array->[1]) > 0) {
		if (defined($kegghash->{$array->[0]}) && defined($compoundhash->{$kegghash->{$array->[0]}})) {
			$compoundhash->{$kegghash->{$array->[0]}}->{search_inchi} = $array->[1];
		}
	}
}
close($fd);
$filename = $Bin."/../Structures/MetaCyc_Charged_InChI.txt";
open(my $fe, "<", $filename);
while (my $line = <$fe>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	if (defined($array->[1]) && length($array->[1]) > 0) {
		if (defined($metacychash->{$array->[0]}) && defined($compoundhash->{$metacychash->{$array->[0]}})) {
			$compoundhash->{$metacychash->{$array->[0]}}->{structure} = $array->[1];
		}
	}
}
close($fe);
$filename = $Bin."/../Structures/MetaCyc_Charged_MolAnalysis.tbl";
open(my $ff, "<", $filename);
while (my $line = <$ff>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	if (defined($array->[1]) && length($array->[1]) > 0) {
		if (defined($metacychash->{$array->[0]}) && defined($compoundhash->{$metacychash->{$array->[0]}})) {
			$compoundhash->{$metacychash->{$array->[0]}}->{groups} = [split(/\|/,$array->[2])];
			$compoundhash->{$metacychash->{$array->[0]}}->{charge} = $array->[3];
			$compoundhash->{$metacychash->{$array->[0]}}->{formula} = $array->[4];
			$compoundhash->{$metacychash->{$array->[0]}}->{mass} = $array->[6];
			$compoundhash->{$metacychash->{$array->[0]}}->{deltag} = $array->[7];
			$compoundhash->{$metacychash->{$array->[0]}}->{deltagerr} = $array->[8];
		}
	}
}
close($ff);

$filename = $Bin."/../Pathways/HopeScenarios.txt";
open(my $fi, "<", $filename);
while (my $line = <$fi>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	if (defined($array->[1]) && length($array->[1]) > 0) {
		if (defined($kegghash->{$array->[1]}) && defined($reactionhash->{$kegghash->{$array->[1]}})) {
			$array->[0] =~ s/:reactions$//;
			push(@{$reactionhash->{$kegghash->{$array->[1]}}->{hope_pathways}},$array->[0]);
		}
	}
}
close($fi);
$filename = $Bin."/../Pathways/plantdefault.pathways.tsv";
open(my $fj, "<", $filename);
$heading = <$fj>;
chomp($heading);
$array = [split(/\t/,$heading)];
$headinghash = {};
for (my $i=0; $i <@{$array}; $i++) {
	$headinghash->{$array->[$i]} = $i;
}
while (my $line = <$fj>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	if (defined($array->[0]) && length($array->[0]) > 0) {
		if (defined($reactionhash->{$array->[0]})) {
			if (defined($headinghash->{KEGG}) && $array->[$headinghash->{KEGG}] ne "null") {
				$reactionhash->{$array->[0]}->{kegg_pathways} = [split(/\|/,$array->[$headinghash->{KEGG}])];
			}
			if (defined($headinghash->{AraCyc}) && $array->[$headinghash->{AraCyc}] ne "null") {
				$reactionhash->{$array->[0]}->{aracyc_pathways} = [split(/\|/,$array->[$headinghash->{AraCyc}])];
			}
			if (defined($headinghash->{EcoCyc}) && $array->[$headinghash->{EcoCyc}] ne "null") {
				$reactionhash->{$array->[0]}->{ecocyc_pathways} = [split(/\|/,$array->[$headinghash->{EcoCyc}])];
			}
			if (defined($headinghash->{MetaCyc}) && $array->[$headinghash->{MetaCyc}] ne "null") {
				$reactionhash->{$array->[0]}->{metacyc_pathways} = [split(/\|/,$array->[$headinghash->{MetaCyc}])];
			}
			if (defined($headinghash->{PlantCyc}) && $array->[$headinghash->{PlantCyc}] ne "null") {
				$reactionhash->{$array->[0]}->{plantcyc_pathways} = [split(/\|/,$array->[$headinghash->{PlantCyc}])];
			}
		}
	}
}
close($fj);

#Reading roles
my $rolehash = {};
my $roleidhash = {};
my $newroleidhash = {};
$filename = $Bin."/../Templates/Roles.tsv";
open(my $fk, "<", $filename);
$heading = <$fk>;
chomp($heading);
$array = [split(/\t/,$heading)];
$headinghash = {};
for (my $i=0; $i <@{$array}; $i++) {
	$headinghash->{$array->[$i]} = $i;
}
while (my $line = <$fk>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	if (defined($array->[0]) && length($array->[0]) > 0) {		
		$rolehash->{$array->[0]} = {
			id => $array->[0],
			name => $array->[$headinghash->{name}],
			source => $array->[$headinghash->{source}],
			complexes => [],
			reactions => [],
			subsystems => [],	
		};
		$rolehash->{$array->[0]}->{searchname} = lc($array->[$headinghash->{name}]);
		$rolehash->{$array->[0]}->{searchname} =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
		$rolehash->{$array->[0]}->{searchname} =~ s/\s//g;
		$rolehash->{$array->[0]}->{searchname} =~ s/\#.*$//g;
		$rolehash->{$array->[0]}->{searchname} =~ s/\(ec\)//g;
		$rolehash->{$array->[0]}->{id} =~ s/ftr0*/fr./;
		$newroleidhash->{$rolehash->{$array->[0]}->{id}} = $rolehash->{$array->[0]};
		if (defined($headinghash->{features}) && $array->[$headinghash->{features}] ne "null" && length($array->[$headinghash->{features}]) > 0) {
			$rolehash->{$array->[0]}->{features} = [split(/;/,$array->[$headinghash->{features}])];
		}
		if (defined($headinghash->{aliases}) && $array->[$headinghash->{aliases}] ne "null" && length($array->[$headinghash->{aliases}]) > 0) {
			$rolehash->{$array->[0]}->{aliases} = [split(/;/,$array->[$headinghash->{aliases}])];
			for (my $i=0; $i < @{$rolehash->{$array->[0]}->{aliases}}; $i++) {
				my $temparray = [split(/:/,$rolehash->{$array->[0]}->{aliases}->[$i])];
				$rolehash->{$array->[0]}->{aliases}->[$i] = pop(@{$temparray});
				if ($rolehash->{$array->[0]}->{aliases}->[$i] =~ m/^msfr\./) {
					$roleidhash->{$rolehash->{$array->[0]}->{aliases}->[$i]} = $rolehash->{$array->[0]};
				}
			}
		}
	}
}
close($fk);

#Reading subsystems
my $sshash = {};
$filename = $Bin."/../mapping.json";
open(my $fm, "<", $filename);
my $mapping = <$fm>;
close($fm);
$mapping = decode_json $mapping;
my $sss = $mapping->{data}->{subsystems};
for (my $i=0; $i < @{$sss}; $i++) {
	$sss->[$i]->{id} =~ s/msrs/ss/;
	$sshash->{$sss->[$i]->{id}} = {
		id => $sss->[$i]->{id},
		type => $sss->[$i]->{type},
		class => $sss->[$i]->{class},
		name => $sss->[$i]->{name},
		subclass => $sss->[$i]->{subclass},
		roles => [],
		reactions => []
	};
	my $array = $sss->[$i]->{role_refs};
	for (my $j=0; $j < @{$array}; $j++) {
		my $subarray = [split(/\//,$array->[$j])];
		my $roleid = pop(@{$subarray});
		if (defined($roleidhash->{$roleid})) {
			push(@{$sshash->{$sss->[$i]->{id}}->{roles}},$roleidhash->{$roleid}->{id}.";".$roleidhash->{$roleid}->{name}.";".$roleidhash->{$roleid}->{searchname});
			push(@{$roleidhash->{$roleid}->{subsystems}},$sshash->{$sss->[$i]->{id}}->{id}.";".$sshash->{$sss->[$i]->{id}}->{name}.";".$sshash->{$sss->[$i]->{id}}->{class}.";".$sshash->{$sss->[$i]->{id}}->{subclass});
		}
	}
}

#Reading complexes
my $cpxhash = {};
$filename = $Bin."/../Templates/Complexes.tsv";
open(my $fl, "<", $filename);
$heading = <$fl>;
chomp($heading);
$array = [split(/\t/,$heading)];
$headinghash = {};
for (my $i=0; $i <@{$array}; $i++) {
	$headinghash->{$array->[$i]} = $i;
}
my $count = 1;
while (my $line = <$fl>) {
	chomp($line);
	$array = [split(/\t/,$line)];
	if (defined($array->[0]) && length($array->[0]) > 0) {
		$cpxhash->{$array->[0]} = {
			id => "cpx.".$count,
			source => $array->[$headinghash->{source}],
			name => "cpx.".$count,
			confidence => $array->[$headinghash->{confidence}],
			roles => [],
			reactions => [],
			aliases => [$array->[0],$array->[$headinghash->{name}]]
		};
		$count++;
		if (defined($headinghash->{reference}) && $array->[$headinghash->{reference}] ne "null" && length($array->[$headinghash->{reference}]) > 0) {
			$cpxhash->{$array->[0]}->{reference} = $array->[$headinghash->{reference}];
		}
		if (defined($headinghash->{roles}) && $array->[$headinghash->{roles}] ne "null" && length($array->[$headinghash->{roles}]) > 0) {
			my $rolearray = [split(/\|/,$array->[$headinghash->{roles}])];
			for (my $i=0; $i < @{$rolearray}; $i++) {
				my $temparray = [split(/;/,$rolearray->[$i])];
				if (defined($rolehash->{$temparray->[0]})) {
					if ($cpxhash->{$array->[0]}->{name} eq $cpxhash->{$array->[0]}->{id}) {
						$cpxhash->{$array->[0]}->{name} = $rolehash->{$temparray->[0]}->{name};
					}
					push(@{$rolehash->{$temparray->[0]}->{complexes}},$cpxhash->{$array->[0]}->{id}.";".$cpxhash->{$array->[0]}->{name});
					push(@{$cpxhash->{$array->[0]}->{roles}},$rolehash->{$temparray->[0]}->{id}.";".$rolehash->{$temparray->[0]}->{name}.";".$rolehash->{$temparray->[0]}->{searchname});
				}
			}
		}
	}
}
close($fl);

#Functional roles of reactions
my $templatelist = ["GramPositive","GramNegative"];
my $rxnrolehash;
my $rxnsshash;
for (my $i=0;$i < @{$templatelist}; $i++) {
	$filename = $Bin."/../Templates/".$templatelist->[$i]."/Reactions.tsv";
	open(my $fl, "<", $filename);
	$heading = <$fl>;
	chomp($heading);
	$array = [split(/\t/,$heading)];
	$headinghash = {};
	for (my $i=0; $i <@{$array}; $i++) {
		$headinghash->{$array->[$i]} = $i;
	}
	while (my $line = <$fl>) {
		chomp($line);
		$array = [split(/\t/,$line)];
		if (defined($array->[0]) && length($array->[0]) > 0) {
			if (defined($reactionhash->{$array->[0]})) {
				if (defined($headinghash->{complexes}) && $array->[$headinghash->{complexes}] ne "null" && length($array->[$headinghash->{complexes}]) > 0) {
					my $cpxarray = [split(/\|/,$array->[$headinghash->{complexes}])];
					for (my $j=0; $j < @{$cpxarray}; $j++) {
						if (defined($cpxhash->{$cpxarray->[$j]})) {
							my $roles = $cpxhash->{$cpxarray->[$j]}->{roles};
							push(@{$cpxhash->{$cpxarray->[$j]}->{reactions}},$array->[0].";".$reactionhash->{$array->[0]}->{definition});
							for (my $k=0; $k < @{$roles}; $k++) {
								my $rolearray = [split(/;/,$roles->[$k])];
								if (defined($newroleidhash->{$rolearray->[0]})) {
									if (!defined($rxnrolehash->{$array->[0]}->{$rolearray->[0]})) {
										push(@{$newroleidhash->{$rolearray->[0]}->{reactions}},$array->[0].";".$reactionhash->{$array->[0]}->{definition});
										$rxnrolehash->{$array->[0]}->{$rolearray->[0]} = 1;
										my $sss = $newroleidhash->{$rolearray->[0]}->{subsystems};
										push(@{$reactionhash->{$array->[0]}->{roles}},$newroleidhash->{$rolearray->[0]}->{id}.";".$newroleidhash->{$rolearray->[0]}->{name}.";".$newroleidhash->{$rolearray->[0]}->{searchname});
										for (my $m=0; $m < @{$sss}; $m++) {
											my $ssarray = [split(/;/,$sss->[$m])];
											if (defined($sshash->{$ssarray->[0]})) {
												if (!defined($rxnsshash->{$array->[0]}->{$ssarray->[0]})) {
													$rxnsshash->{$array->[0]}->{$ssarray->[0]} = 1;
													push(@{$reactionhash->{$array->[0]}->{subsystems}},$sshash->{$ssarray->[0]}->{id}.";".$sshash->{$ssarray->[0]}->{name}.";".$sshash->{$ssarray->[0]}->{class}.";".$sshash->{$ssarray->[0]}->{subclass});
													push(@{$sshash->{$ssarray->[0]}->{reactions}},$array->[0].";".$reactionhash->{$array->[0]}->{definition});
												}
											}
										}
									}
								}
							}
							push(@{$reactionhash->{$array->[0]}->{complexes}},$cpxhash->{$cpxarray->[$j]}->{id}.";".$cpxhash->{$cpxarray->[$j]}->{name});
						}
					}
					push(@{$reactionhash->{$array->[0]}->{templates}},$templatelist->[$i].";".$array->[$headinghash->{type}].";".$array->[$headinghash->{compartment}].";".$array->[$headinghash->{direction}].";".$array->[$headinghash->{gfdir}].";".$array->[$headinghash->{base_cost}].";".$array->[$headinghash->{forward_cost}].";".$array->[$headinghash->{reverse_cost}]);
				}
			}
		}
	}
}

my $JSON = JSON::XS->new->ascii->pretty;

my $jsonobject = [];
foreach my $key (keys(%{$reactionhash})) {
	foreach my $field (keys(%{$reactionhash->{$key}})) {
		if (ref($reactionhash->{$key}->{$field}) eq "ARRAY") {
			@{$reactionhash->{$key}->{$field}} = uniq(@{$reactionhash->{$key}->{$field}});
		}
	}
	push(@{$jsonobject},$reactionhash->{$key});
}
$filename = $Bin."/../SOLRDump/Reactions.json";
open(my $fn, ">", $filename);
print $fn $JSON->encode($jsonobject);
close($fn);

$jsonobject = [];
foreach my $key (keys(%{$compoundhash})) {
	foreach my $field (keys(%{$compoundhash->{$key}})) {
		if (ref($compoundhash->{$key}->{$field}) eq "ARRAY") {
			@{$compoundhash->{$key}->{$field}} = uniq(@{$compoundhash->{$key}->{$field}});
		}
	}
	push(@{$jsonobject},$compoundhash->{$key});
}
$filename = $Bin."/../SOLRDump/Compounds.json";
open(my $fo, ">", $filename);
print $fo $JSON->encode($jsonobject);
close($fo);

$jsonobject = [];
foreach my $key (keys(%{$rolehash})) {
	foreach my $field (keys(%{$rolehash->{$key}})) {
		if (ref($rolehash->{$key}->{$field}) eq "ARRAY") {
			@{$rolehash->{$key}->{$field}} = uniq(@{$rolehash->{$key}->{$field}});
		}
	}
	push(@{$jsonobject},$rolehash->{$key});
}
$filename = $Bin."/../SOLRDump/Roles.json";
open(my $fp, ">", $filename);
print $fp $JSON->encode($jsonobject);
close($fp);

$jsonobject = [];
foreach my $key (keys(%{$cpxhash})) {
	foreach my $field (keys(%{$cpxhash->{$key}})) {
		if (ref($cpxhash->{$key}->{$field}) eq "ARRAY") {
			@{$cpxhash->{$key}->{$field}} = uniq(@{$cpxhash->{$key}->{$field}});
		}
	}
	push(@{$jsonobject},$cpxhash->{$key});
}
$filename = $Bin."/../SOLRDump/Complexes.json";
open(my $fq, ">", $filename);
print $fq $JSON->encode($jsonobject);
close($fq);

$jsonobject = [];
foreach my $key (keys(%{$sshash})) {
	foreach my $field (keys(%{$sshash->{$key}})) {
		if (ref($sshash->{$key}->{$field}) eq "ARRAY") {
			@{$sshash->{$key}->{$field}} = uniq(@{$sshash->{$key}->{$field}});
		}
	}
	push(@{$jsonobject},$sshash->{$key});
}
$filename = $Bin."/../SOLRDump/Subsystems.json";
open(my $fr, ">", $filename);
print $fr $JSON->encode($jsonobject);
close($fr);

print "Primary class\tSub class\tSubsystem\tRole\n";
for (my $i=0; $i < @{$jsonobject}; $i++) {
	if ($jsonobject->[$i]->{class} eq "Cofactors, Vitamins, Prosthetic Groups, Pigments") {
		for (my $j=0; $j < @{$jsonobject->[$i]->{roles}}; $j++) {
			print $jsonobject->[$i]->{class}."\t".$jsonobject->[$i]->{subclass}."\t".$jsonobject->[$i]->{name}."\t".$jsonobject->[$i]->{roles}->[$j]."\n";
		}
	}
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

1;
