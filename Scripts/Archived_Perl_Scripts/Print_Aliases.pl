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

my %BiGG = (iAF1260=>1,iAF692=>1,iIN800=>1,iJR904=>1,
	    iMA945=>1,iMEO21=>1,iMM904=>1,iRR1083=>1,
	    iSB619=>1,iSO783=>1);
	    
my %Global_Aliases=();
my %Cpd_Aliases = %{$Default_JSON->{"compound_aliases"}};
foreach my $cpd (keys %Cpd_Aliases){
    foreach my $alias (keys %{$Cpd_Aliases{$cpd}}){
	foreach my $entry (@{$Cpd_Aliases{$cpd}{$alias}}){
	    $entry =~ s/^\s+//;
	    $entry =~ s/\s+$//;

	    $Global_Aliases{"Compounds"}{$alias}{$entry}{default}{$cpd}=1;
	    if(exists($BiGG{$alias})){
		$Global_Aliases{"Compounds"}{"BiGG"}{$entry}{default}{$cpd}=1;
	    }
	}
    }
}

my %Rxn_Aliases = %{$Default_JSON->{"reaction_aliases"}};
foreach my $rxn (keys %Rxn_Aliases){
    foreach my $alias (keys %{$Rxn_Aliases{$rxn}}){
	foreach my $entry (@{$Rxn_Aliases{$rxn}{$alias}}){
	    $entry =~ s/^\s+//;
	    $entry =~ s/\s+$//;

	    if($alias eq "Enzyme Class"){
		foreach my $sub_entry (split(/[\s\/]+/,$entry)){
		    next if $sub_entry =~ /determined/;
		    $sub_entry =~ s/^EC-//;
		    $sub_entry =~ s/[\),]+//;

		    $Global_Aliases{"Reactions"}{$alias}{$sub_entry}{default}{$rxn}=1;
		    if(exists($BiGG{$alias})){
			$Global_Aliases{"Reactions"}{"BiGG"}{$sub_entry}{default}{$rxn}=1;
		    }
		}
	    }else{
		$Global_Aliases{"Reactions"}{$alias}{$entry}{default}{$rxn}=1;
		if(exists($BiGG{$alias})){
		    $Global_Aliases{"Reactions"}{"BiGG"}{$entry}{default}{$rxn}=1;
		}
	    }
	}
    }
}

%Cpd_Aliases = %{$PlantDefault_JSON->{"compound_aliases"}};
foreach my $cpd (keys %Cpd_Aliases){
    foreach my $alias ( grep { $_ ne "ModelSEED" } 
			keys %{$Cpd_Aliases{$cpd}}){
	foreach my $entry (@{$Cpd_Aliases{$cpd}{$alias}}){
	    $entry =~ s/^\s+//;
	    $entry =~ s/\s+$//;

	    $Global_Aliases{"Compounds"}{$alias}{$entry}{plantdefault}{$cpd}=1;
	    if(exists($BiGG{$alias})){
		$Global_Aliases{"Compounds"}{"BiGG"}{$entry}{plantdefault}{$cpd}=1;
	    }
	}
    }
}

%Rxn_Aliases = %{$PlantDefault_JSON->{"reaction_aliases"}};
foreach my $rxn (keys %Rxn_Aliases){
    foreach my $alias ( grep { $_ ne "ModelSEED" } 
			keys %{$Rxn_Aliases{$rxn}}){
	foreach my $entry (@{$Rxn_Aliases{$rxn}{$alias}}){
	    $entry =~ s/^\s+//;
	    $entry =~ s/\s+$//;

	    if($alias eq "Enzyme Class"){
		foreach my $sub_entry (split(/[\s\/]+/,$entry)){
		    next if $sub_entry =~ /determined/;
		    $sub_entry =~ s/^EC-//;
		    $sub_entry =~ s/[\),]+//;

		    $Global_Aliases{"Reactions"}{$alias}{$sub_entry}{plantdefault}{$rxn}=1;
		    if(exists($BiGG{$alias})){
			$Global_Aliases{"Reactions"}{"BiGG"}{$sub_entry}{plantdefault}{$rxn}=1;
		    }
		}
	    }else{
		$Global_Aliases{"Reactions"}{$alias}{$entry}{plantdefault}{$rxn}=1;
		if(exists($BiGG{$alias})){
		    $Global_Aliases{"Reactions"}{"BiGG"}{$entry}{plantdefault}{$rxn}=1;
		}
	    }
	}
    }
}

foreach my $entity (["Compounds","cpd"],["Reactions","rxn"]){
    my $file = $entity->[0]."_Aliases.tsv";
    open(MAIN, "> ".$file);
    print MAIN "MS ID\tOld MS ID\tExternal ID\tSource\n";

    open(NAMES, "> Names_".$file);
    print NAMES "MS ID\tOld MS ID\tExternal ID\tSource\n";

    open(MODELS, "> Models_".$file);
    print MODELS "MS ID\tOld MS ID\tExternal ID\tSource\n";

    open(BIOCYC, "> BioCyc_".$file);
    print BIOCYC "MS ID\tOld MS ID\tExternal ID\tSource\n";

    open(EC, "> Enzyme_Class_".$file);
    print EC "MS ID\tOld MS ID\tExternal ID\tSource\n";

    foreach my $aliasSet (sort keys %{$Global_Aliases{$entity->[0]}}){
	if($aliasSet =~ /KEGG|MetaCyc|BiGG|PlantCyc/){
	    foreach my $alias (sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}}){
		print MAIN exists($Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{plantdefault}) ? join("|",sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{plantdefault}}) : "";
		print MAIN "\t";
		print MAIN exists($Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{default}) ? join("|",sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{default}}) : "";
		print MAIN "\t";
		print MAIN $alias."\t".$aliasSet."\n";
	    }
	}elsif($aliasSet =~ /Cyc$/ && $aliasSet !~ /BioCyc/){
	    foreach my $alias (sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}}){
		print BIOCYC exists($Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{plantdefault}) ? join("|",sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{plantdefault}}) : "";
		print BIOCYC "\t";
		print BIOCYC exists($Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{default}) ? join("|",sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{default}}) : "";
		print BIOCYC "\t";
		print BIOCYC $alias."\t".$aliasSet."\n";
	    }
	}elsif($aliasSet eq "Enzyme Class"){
	    foreach my $alias (sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}}){
		print EC exists($Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{plantdefault}) ? join("|",sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{plantdefault}}) : "";
		print EC "\t";
		print EC exists($Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{default}) ? join("|",sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{default}}) : "";
		print EC "\t";
		print EC $alias."\t".$aliasSet."\n";
	    }
	}elsif($aliasSet =~ /name/){
	    foreach my $alias (sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}}){
		print NAMES exists($Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{plantdefault}) ? join("|",sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{plantdefault}}) : "";
		print NAMES "\t";
		print NAMES exists($Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{default}) ? join("|",sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{default}}) : "";
		print NAMES "\t";
		print NAMES $alias."\t".$aliasSet."\n";
	    }
	}elsif($aliasSet !~ /KBase|ModelSEED|MicrobesOnlineFitness|NewMediaTransporters|obsolete|BioCyc/){
	    foreach my $alias (sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}}){
		print MODELS exists($Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{plantdefault}) ? join("|",sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{plantdefault}}) : "";
		print MODELS "\t";
		print MODELS exists($Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{default}) ? join("|",sort keys %{$Global_Aliases{$entity->[0]}{$aliasSet}{$alias}{default}}) : "";
		print MODELS "\t";
		print MODELS $alias."\t".$aliasSet."\n";
	    }
	}
    }
    close(MAIN);
    close(BIOCYC);
    close(NAMES);
    close(MODELS);
    close(EC);
}
