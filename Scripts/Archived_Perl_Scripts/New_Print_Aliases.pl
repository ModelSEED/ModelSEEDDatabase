#!/usr/bin/env perl
use warnings;
use strict;
use JSON;

open(FH, "< default.obj");
my $Default_JSON;
while(<FH>){
    $Default_JSON.=$_;
}
close(FH);
$Default_JSON=from_json($Default_JSON);

open(FH, "< plantdefault_obs.obj");
my $PlantDefault_JSON;
while(<FH>){
    $PlantDefault_JSON.=$_;
}
close(FH);
$PlantDefault_JSON=from_json($PlantDefault_JSON);

my %Global_Aliases=();
my %Cpd_Aliases = %{$Default_JSON->{"compound_aliases"}};
foreach my $cpd (keys %Cpd_Aliases){
    foreach my $alias (keys %{$Cpd_Aliases{$cpd}}){
	foreach my $entry (@{$Cpd_Aliases{$cpd}{$alias}}){
	    $Global_Aliases{$alias}{$entry}{default}{$cpd}=1;
	}
    }
}

my %Rxn_Aliases = %{$Default_JSON->{"reaction_aliases"}};
foreach my $rxn (keys %Rxn_Aliases){
    foreach my $alias (keys %{$Rxn_Aliases{$rxn}}){
	foreach my $entry (@{$Rxn_Aliases{$rxn}{$alias}}){
	    $Global_Aliases{$alias}{$entry}{default}{$rxn}=1;
	}
    }
}

%Cpd_Aliases = %{$PlantDefault_JSON->{"compound_aliases"}};
foreach my $cpd (keys %Cpd_Aliases){
    foreach my $alias ( grep { $_ ne "ModelSEED" } 
			keys %{$Cpd_Aliases{$cpd}}){
	foreach my $entry (@{$Cpd_Aliases{$cpd}{$alias}}){
	    $Global_Aliases{$alias}{$entry}{plantdefault}{$cpd}=1;
	}
    }
}

foreach my $cpd (keys %Cpd_Aliases){
    if(exists($Cpd_Aliases{$cpd}{"ModelSEED"})){
	foreach my $link (@{$Cpd_Aliases{$cpd}{"ModelSEED"}}){
	    
	    foreach my $alias (keys %{$Cpd_Aliases{$link}}){
		foreach my $entry (@{$Cpd_Aliases{$link}{$alias}}){
		    $Global_Aliases{$alias}{$entry}{plantdefault}{$cpd}=1;
		}
	    }

	    foreach my $alias (keys %{$Cpd_Aliases{$cpd}}){
		foreach my $entry (@{$Cpd_Aliases{$cpd}{$alias}}){
		    $Global_Aliases{$alias}{$entry}{plantdefault}{$link}=1;
		}
	    }
	}
    }
}

%Rxn_Aliases = %{$PlantDefault_JSON->{"reaction_aliases"}};
foreach my $rxn (keys %Rxn_Aliases){
    foreach my $alias (keys %{$Rxn_Aliases{$rxn}}){
	foreach my $entry (@{$Rxn_Aliases{$rxn}{$alias}}){
	    $Global_Aliases{$alias}{$entry}{plantdefault}{$rxn}=1;
	}
    }
}

foreach my $aliasSet (sort keys %Global_Aliases){
    my $file = $aliasSet;
    $file = join("_",split(/\s/,$aliasSet)) if $aliasSet eq "Enzyme Class";
    $file = "../Aliases/".$file.".aliases";

    open(OUT, "> ".$file);
    print OUT $aliasSet."\tdefault\tplantdefault\n";
    foreach my $alias (sort keys %{$Global_Aliases{$aliasSet}}){
	print OUT $alias."\t";
	print OUT exists($Global_Aliases{$aliasSet}{$alias}{default}) ? join("|",sort keys %{$Global_Aliases{$aliasSet}{$alias}{default}}) : "";
	print OUT "\t";
	print OUT exists($Global_Aliases{$aliasSet}{$alias}{plantdefault}) ? join("|",sort keys %{$Global_Aliases{$aliasSet}{$alias}{plantdefault}}) : ""; 
	print OUT "\n";
    }
    close(OUT);
}
