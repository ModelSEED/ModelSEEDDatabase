#!/usr/bin/env perl
#use PlantSEED::Formulas;
use warnings;
use strict;
my @temp=();

my $Database_Root = "/homes/seaver/Biochemistry_Mirrors/ftp.bioinformatics.jp/kegg/ligand/";
my $Structure_Root = "/homes/seaver/Projects/ModelSEEDDatabase/Biochemistry/Structures/KEGG/";
my $Output_Root = "/homes/seaver/Projects/ModelSEEDDatabase/Biochemistry/Aliases/Provenance/Primary_Databases/";

my %Structures=();
foreach my $type ("InChI","SMILE"){
    open(FH, "< ".$Structure_Root.$type."_ChargedStrings.txt");
    while(<FH>){
	chomp;
	@temp=split(/\t/,$_);
	$Structures{$temp[0]}{$type}=$temp[1];
    }

    open(FH, "< ".$Structure_Root.$type."_OriginalStrings.txt");
    while(<FH>){
	chomp;
	@temp=split(/\t/,$_);
	next if exists($Structures{$temp[0]}) && exists($Structures{$temp[0]}{$type});
	$Structures{$temp[0]}{$type}=$temp[1];
    }
}

open(FH, "< ".$Database_Root."compound/compound");
my ($field,$data)="";
my $ENTRY="";
my @NAMES=();
my $FORMULA="null";
my $MASS=10000000;
my $CHARGE=10000000;
my $Compounds={};
while(<FH>){
    chomp;
    @temp=split /\s+/;
    if($temp[0] ne ""){
	($field,$data)=($temp[0],join(" ",@temp[1..$#temp]));
    }else{
	$data=join(" ",@temp[1..$#temp]);
    }

    if($field eq "ENTRY"){
	$ENTRY=$temp[1];
    }

    if($field eq "NAME"){
	$data =~ s/;$//;
	push(@NAMES,$data);
    }

    if($field eq "FORMULA"){
#	unless($data =~ /\)n/ || ($data =~ /n/ && $data !~ /[MZSIR]n/)){
	    $FORMULA=$data;
#	}
    }

    if($field eq "EXACT_MASS"){
	$MASS=$data;
    }

    if($temp[0] eq "///"){

	#Adjust charge and mass for Generic "R" groups
	if($FORMULA =~ /R/ && $FORMULA ne "Rn"){
	    $CHARGE = 0;
	    $MASS = 0;
	}

	$Compounds->{$ENTRY}={Compounds=>$ENTRY,NAMES=>join("|",@NAMES),FORMULA=>$FORMULA,CHARGE=>$CHARGE,MASS=>$MASS};

	#Modifictions for electron/photon
	if($ENTRY eq "C05359"){
	    $Compounds->{$ENTRY}{MASS}=0;
	    $Compounds->{$ENTRY}{CHARGE}=-1;
	}
	
	if($ENTRY eq "C00205"){
	    $Compounds->{$ENTRY}{MASS}=0;
	    $Compounds->{$ENTRY}{CHARGE}=0;
	}
	
	($ENTRY,$FORMULA,$MASS,$CHARGE)=("","null",10000000,10000000);
	undef(@NAMES);
    }
}
close(FH);

open(FH, "< ".$Database_Root."glycan/glycan");
($field,$data)=("","");
($ENTRY,$FORMULA)=("","null");
my $Glycans={};
undef(@NAMES);
undef($Glycans);
while(<FH>){
    chomp;
    @temp=split /\s+/;
    if($temp[0] ne ""){
	($field,$data)=($temp[0],join(" ",@temp[1..$#temp]));
    }else{
	$data=join(" ",@temp[1..$#temp]);
    }
    
    if($field eq "ENTRY"){
	$ENTRY=$temp[1];
    }
    
    if($field eq "NAME"){
	$data =~ s/;$//;
	push(@NAMES,$data);
    }

    if($field eq "COMPOSITION" && scalar(@NAMES)==0){
	$data =~ s/\s//g;
	push(@NAMES,$data);
    }

    if($field eq "FORMULA"){
	$FORMULA=$data;
    }

    if($temp[0] eq "///"){
	$Glycans->{$ENTRY}={Glycans=>$ENTRY,NAMES=>join("|",@NAMES),FORMULA=>$FORMULA,CHARGE=>10000000,MASS=>10000000};
	($ENTRY,$FORMULA)=("","null");
	undef(@NAMES);
    }
}
close(FH);

my $Pathways={};
open(FH, "< ".$Database_Root."../pathway/pathway");
($field,$data,$ENTRY)=("","","");
my ($CLASS,$NAME)=("","","");
while(<FH>){
    chomp;
    @temp=split /\s+/;
    if($temp[0] ne ""){
	($field,$data)=($temp[0],join(" ",@temp[1..$#temp]));
    }else{
	$data=join(" ",@temp[1..$#temp]);
    }

    if($field eq "ENTRY"){
	$ENTRY=$temp[1];
    }
    if($field eq "NAME"){
	$NAME=$data;
    }
    if($field eq "CLASS"){
	$CLASS=$data;
    }
    if($temp[0] eq "///"){
	next if $ENTRY !~ /^map\d+/;
	$ENTRY =~ s/map/rn/;
	$Pathways->{$ENTRY}->{"C"}=$CLASS;
	$Pathways->{$ENTRY}->{"N"}=$NAME;
	($ENTRY,$CLASS,$NAME)=("","","");
    }
}
close(FH);

open(FH, "< ".$Database_Root."reaction/reaction");
undef(@temp);
($field,$data,$ENTRY)=("","","");
undef(@NAMES);
my $EQUATION="";
my @ECs=();

open(OUT, "> KEGG_Reaction_Pathway_Mapping.txt");
my @Pathways=();

my $Reactions={};
my $count=0;
while(<FH>){
    chomp;
    @temp=split /\s+/;
    if($temp[0] ne ""){
	($field,$data)=($temp[0],join(" ",@temp[1..$#temp]));
    }else{
	$data=join(" ",@temp[1..$#temp]);
    }

    if($field eq "ENTRY"){
	$ENTRY=$temp[1];
    }

    if($field eq "NAME"){
	$data =~ s/;$//;
	push(@NAMES,$data);
    }

    if($field eq "EQUATION"){
	$EQUATION=$data;
    }

    if($field eq "ENZYME"){
	push(@ECs, split(/\s/,$data));
    }

    if($field eq "PATHWAY"){
	next if $data =~ /rn01100/; #skip generic "Metabolic pathways"
	next if $data =~ /rn01110/; #skip generic "Biosynthesis of secondary metabolites"
	next if $data =~ /rn01120/; #skip generic "Microbial metabolism in diverse environments"
	push(@Pathways,(split(/\s+/,$data))[0]);
    }

    if($field eq "///"){
	$Reactions->{$count}={"KEGG"=>$ENTRY,"NAMES"=>"".join("|",@NAMES),"EQUATION"=>convert_KEGG_equation($EQUATION),"ECs"=>"".join("|",@ECs)};
	$Reactions->{$count}{'ECs'}="null" if !$Reactions->{$count}{'ECs'};

	while($EQUATION =~ /([CG]\d{5})/g){
	    if(!exists($Compounds->{$1})){
		if(exists($Glycans->{$1})){
		    $Compounds->{$1}=$Glycans->{$1};
		}else{
		    print "Missing $1 in $ENTRY\n";
		}
	    }
	}

	foreach my $p (@Pathways){
	    print OUT $ENTRY,"\t";
	    print OUT join("|",@ECs),"\t";
	    print OUT $p,"\t";
	    print OUT $Pathways->{$p}->{"N"},"\t";
	    print OUT $Pathways->{$p}->{"C"},"\n";
	}

	($ENTRY,$EQUATION)=("","");
	undef(@NAMES);
	undef(@ECs);
	undef(@Pathways);
    }
    $count++;
}
close(FH);
close(OUT);

open(OUT, "> ".$Output_Root."KEGG_Compounds.tbl");
my @Headers=("NAMES","FORMULA","CHARGE","MASS");
print OUT "ID\t".join("\t",@Headers)."\tInChI\tSMILE\n";
foreach my $id (sort keys %$Compounds){
    print OUT $id,"\t";
    foreach my $h (@Headers){
	print OUT $Compounds->{$id}{$h};
	print OUT "\t"; # unless $h eq $Headers[$#Headers];
    }
    foreach my $type ("InChI","SMILE"){
	if(exists($Structures{$id}) && exists($Structures{$id}{$type})){
	    print OUT $Structures{$id}{$type};
	}else{
	    print OUT "-";
	}
	print OUT "\t";
    }
    print OUT "\n";
}
close(OUT);

open(OUT, "> ".$Output_Root."KEGG_Glycans.tbl");
print OUT "ID\t".join("\t",@Headers)."\n";
foreach my $id (sort keys %$Glycans){
    print OUT $id,"\t";
    foreach my $h (@Headers){
	print OUT $Glycans->{$id}{$h};
	print OUT "\t" unless $h eq $Headers[$#Headers];
    }
    print OUT "\n";
}
close(OUT);

open(OUT, "> ".$Output_Root."KEGG_Reactions.tbl");
@Headers = ("NAMES","KEGG","EQUATION","ECs");
print OUT "ID\t".join("\t",@Headers)."\n";
foreach my $c (sort { $a <=> $b } keys %$Reactions){
    my @ids=sort grep { substr($_,0,1) eq "R" } split /\|/,$Reactions->{$c}{"KEGG"};

    print OUT $ids[0],"\t";
    foreach my $h (@Headers){
	if(!defined($Reactions->{$c}{$h})){
	    print STDERR "Missing $h\n";
	}
	print OUT $Reactions->{$c}{$h};
	print OUT "\t" unless $h eq $Headers[$#Headers];
    }
    print OUT "\n";
}
close(OUT);

sub convert_KEGG_equation{
    my $eq=shift;

    #Four situations
    #1) No modifier
    #2) +1
    #3) -1
    #4) *
    my $n = "";
    my $m = "";
    if($eq =~ /n\-/){
	$n = 2;

	while($eq =~ /(\d)n/g){
	    my $mod = $n * $1;
	    $eq =~ s/$1n/$mod/;
	}

	$eq =~ s/n/$n/g;
	$eq =~ s/2\-(\d)/2-$1/eg;
	$eq =~ s/2\+(\d)/2+$1/eg;
    }

    if($eq =~ /m\-/){
	$m = 2;

	while($eq =~ /(\d)m/g){
	    my $mod = $m * $1;
	    $eq =~ s/$1m/$mod/;
	}

	$eq =~ s/m/$m/g;
	$eq =~ s/2\-(\d)/2-$1/eg;
	$eq =~ s/2\+(\d)/2+$1/eg;
    }

    if($eq =~ /(n\+)|(\+n)/){
	$n = 1;

	while($eq =~ /(\d)n/g){
	    my $mod = $n * $1;
	    $eq =~ s/$1n/$mod/;
	}

	$eq =~ s/n/$n/g;
	$eq =~ s/1\+(\d)/$1+1/eg;
    }

    if($eq =~ /(m\+)|(\+m)/){
	$m = 1;

	while($eq =~ /(\d)m/g){
	    my $mod = $m * $1;
	    $eq =~ s/$1m/$mod/;
	}

	$eq =~ s/m/$m/g;
	$eq =~ s/1\+(\d)/$1+1/eg;
    }

    if($eq =~ /\d[nm]/){
	$eq =~ s/(\d)[nm]/$1/g;
    }

    $eq =~ s/[nm]/1/g;

    my ($left,$right)=split "<=>",$eq;

    my $new_eq="";
    my @cmpds=split /\+/,$left;
    for(my $i=0;$i<scalar(@cmpds);$i++){
	$cmpds[$i] =~ s/^\s+//;
	$cmpds[$i] =~ s/\s+$//;

	#change position of bracketed stoichs
	$cmpds[$i] =~ s/([CG]\d{5})\((\d)\)/$2 $1/;

	my @stoich=split /\s/,$cmpds[$i];
	if(scalar(@stoich)==2){
            if($stoich[0] !~ /^\(\d+\)$/ && $stoich[0] ne "1"){
                $cmpds[$i]="(".$stoich[0].") ".$stoich[1];
            }elsif($stoich[0] ne "1"){
                $cmpds[$i]=$stoich[0]." ".$stoich[1];
            }else{
                $cmpds[$i]=$stoich[1];
            }
	}
    }

    #Add compartment
    for(my $i=0;$i<scalar(@cmpds);$i++){
	$cmpds[$i].="[0]";
    }

    $new_eq=join(" + ",@cmpds)." <=> ";

    @cmpds=split /\+/,$right;
#    print $right,"R\n";
    for(my $i=0;$i<scalar(@cmpds);$i++){
	$cmpds[$i] =~ s/^\s+//;
	$cmpds[$i] =~ s/\s+$//;

	#change position of bracketed stoichs
	$cmpds[$i] =~ s/([CG]\d{5})\((\d)\)/$2 $1/;

	my @stoich=split /\s/,$cmpds[$i];
	if(scalar(@stoich)==2){
            if($stoich[0] !~ /^\(\d+\)$/ && $stoich[0] ne "1"){
                $cmpds[$i]="(".$stoich[0].") ".$stoich[1];
            }elsif($stoich[0] ne "1"){
                $cmpds[$i]=$stoich[0]." ".$stoich[1];
            }else{
                $cmpds[$i]=$stoich[1];
            }
	}
    }

    #Add compartment
    for(my $i=0;$i<scalar(@cmpds);$i++){
	$cmpds[$i].="[0]";
    }

    $new_eq.=join(" + ",@cmpds);

#    print $new_eq,"\n";
    return $new_eq;
}
