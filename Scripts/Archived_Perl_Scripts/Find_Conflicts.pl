#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my @temp=();
my %Source_Conflicts=();

my $BPath = "../../Biochemistry/";
open(FH, "< ".$BPath."Aliases/Compounds_Aliases.tsv");
while(<FH>){
    chomp;
    my ($pmsid,$dmsid,$sid,$source)=split(/\t/,$_,-1);
    next unless $source =~ /^(KEGG|BiGG)$/;

    $Source_Conflicts{$source}{"compounds"}{$sid}{'p'}{$pmsid}=1;
    $Source_Conflicts{$source}{"compounds"}{$sid}{'d'}{$dmsid}=1;
}
close(FH);

open(FH, "< ".$BPath."Aliases/Reactions_Aliases.tsv");
while(<FH>){
    chomp;
    my ($pmsid,$dmsid,$sid,$source)=split(/\t/,$_,-1);
    next unless $source =~ /^(KEGG|BiGG)$/;

    $Source_Conflicts{$source}{"reactions"}{$sid}{'p'}{$pmsid}=1;
    $Source_Conflicts{$source}{"reactions"}{$sid}{'d'}{$dmsid}=1;
}
close(FH);

foreach my $source (sort keys %Source_Conflicts){
    foreach my $type (sort keys %{$Source_Conflicts{$source}}){
	foreach my $id (sort keys %{$Source_Conflicts{$source}{$type}}){
	    my $pstring = join("|", sort keys %{$Source_Conflicts{$source}{$type}{$id}{'p'}});
	    my $dstring = join("|", sort keys %{$Source_Conflicts{$source}{$type}{$id}{'d'}});
	    if($pstring && $dstring && $pstring ne $dstring){
		print $source,"\t",$type,"\t",$id,"\t",$pstring,"\t",$dstring,"\n";
	    }
	}
    }
}

__END__
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
	    $Global_Aliases{"Reactions"}{$alias}{$entry}{default}{$rxn}=1;
	    if(exists($BiGG{$alias})){
		$Global_Aliases{"Reactions"}{"BiGG"}{$entry}{default}{$rxn}=1;
	    }
	}
    }
}

%Cpd_Aliases = %{$PlantDefault_JSON->{"compound_aliases"}};
foreach my $cpd (keys %Cpd_Aliases){
    foreach my $alias ( grep { $_ ne "ModelSEED" } 
			keys %{$Cpd_Aliases{$cpd}}){
	foreach my $entry (@{$Cpd_Aliases{$cpd}{$alias}}){
	    $Global_Aliases{"Compounds"}{$alias}{$entry}{plantdefault}{$cpd}=1;
	    if(exists($BiGG{$alias})){
		$Global_Aliases{"Compounds"}{"BiGG"}{$entry}{plantdefault}{$cpd}=1;
	    }
	}
    }
}

#foreach my $cpd (keys %Cpd_Aliases){
#    if(exists($Cpd_Aliases{$cpd}{"ModelSEED"})){
#	foreach my $link (@{$Cpd_Aliases{$cpd}{"ModelSEED"}}){
	    
#	    foreach my $alias (keys %{$Cpd_Aliases{$link}}){
#		foreach my $entry (@{$Cpd_Aliases{$link}{$alias}}){
#		    $Global_Aliases{"Compounds"}{$alias}{$entry}{plantdefault}{$cpd}=1;
#		    if(exists($BiGG{$alias})){
#			$Global_Aliases{"Compounds"}{"BiGG"}{$entry}{plantdefault}{$cpd}=1;
#		    }
#		}
#	    }

#	    foreach my $alias (keys %{$Cpd_Aliases{$cpd}}){
#		foreach my $entry (@{$Cpd_Aliases{$cpd}{$alias}}){
#		    $Global_Aliases{"Compounds"}{$alias}{$entry}{plantdefault}{$link}=1;
#		    if(exists($BiGG{$alias})){
#			$Global_Aliases{"Compounds"}{"BiGG"}{$entry}{plantdefault}{$link}=1;
#		    }
#		}
#	    }
#	}
#    }
#}

%Rxn_Aliases = %{$PlantDefault_JSON->{"reaction_aliases"}};
foreach my $rxn (keys %Rxn_Aliases){
    foreach my $alias ( grep { $_ ne "ModelSEED" } 
			keys %{$Rxn_Aliases{$rxn}}){
	foreach my $entry (@{$Rxn_Aliases{$rxn}{$alias}}){
	    $Global_Aliases{"Reactions"}{$alias}{$entry}{plantdefault}{$rxn}=1;
	    if(exists($BiGG{$alias})){
		$Global_Aliases{"Reactions"}{"BiGG"}{$entry}{plantdefault}{$rxn}=1;
	    }
	}
    }
}

#foreach my $rxn (keys %Rxn_Aliases){
#    if(exists($Rxn_Aliases{$rxn}{"ModelSEED"})){
#	foreach my $link (@{$Rxn_Aliases{$rxn}{"ModelSEED"}}){
	    
#	    foreach my $alias (keys %{$Rxn_Aliases{$link}}){
#		foreach my $entry (@{$Rxn_Aliases{$link}{$alias}}){
#		    $Global_Aliases{"Reactions"}{$alias}{$entry}{plantdefault}{$rxn}=1;
#		    if(exists($BiGG{$alias})){
#			$Global_Aliases{"Reactions"}{"BiGG"}{$entry}{plantdefault}{$rxn}=1;
#		    }
#		}
#	    }

#	    foreach my $alias (keys %{$Rxn_Aliases{$rxn}}){
#		foreach my $entry (@{$Rxn_Aliases{$rxn}{$alias}}){
#		    $Global_Aliases{"Reactions"}{$alias}{$entry}{plantdefault}{$link}=1;
#		    if(exists($BiGG{$alias})){
#			$Global_Aliases{"Reactions"}{"BiGG"}{$entry}{plantdefault}{$link}=1;
#		    }
#		}
#	    }
#	}
#    }
#}


foreach my $entity (["Compounds","cpd"],["Reactions","rxn"]){
    my $file = $entity->[0]."_Aliases.tsv";
    open(MAIN, "> ".$file);
    print MAIN "MS ID\tOld MS ID\tExternal ID\tSource\n";

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
	}elsif($aliasSet !~ /KBase|ModelSEED|MicrobesOnlineFitness|name|NewMediaTransporters|obsolete|BioCyc/){
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
    close(MODELS);
    close(EC);
}
