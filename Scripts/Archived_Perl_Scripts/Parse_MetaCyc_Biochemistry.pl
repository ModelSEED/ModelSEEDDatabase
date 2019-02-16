#!/usr/bin/env perl
use PlantSEED::Genomes;
use PlantSEED::InChIs;
use PlantSEED::MolAnalyses;
use PlantSEED::Formulas;
use Storable qw(dclone);
use warnings;
use strict;
my @temp=();

sub translate_html_entities{
    my $bit=shift;

    #remove greek letters
    $bit =~ s/\&[Aa]lpha;?/alpha/g;
    $bit =~ s/\&[Bb]eta;?/beta/g;
    $bit =~ s/\&[Gg]amma;/gamma/g;
    $bit =~ s/\&[Dd]elta;/delta/g;
    $bit =~ s/\&[Ee]psilon;?/epsilon/g;
    $bit =~ s/\&[Oo]mega;/omega/g;
    $bit =~ s/\&[Mm]u;/mu/g;
    $bit =~ s/\&[Nn]u;/nu/g;
    $bit =~ s/\&[Kk]appa;/kappa/g;
    $bit =~ s/\&[Cc]hi;/chi/g;
    $bit =~ s/\&[Zz]eta;/zeta/g;
    $bit =~ s/\&[Pp]si;/psi/g;
    $bit =~ s/\&[Pp]i;/pi/g;
    $bit =~ s/\&[Pp]hi;/phi/g;
    $bit =~ s/\&[Tt]au;/tau/g;
    $bit =~ s/\&[Ii]ota;/iota/g;
    $bit =~ s/\&[Tt]heta;/theta/g;
    $bit =~ s/\&[Ss]igma;/sigma/g;
    $bit =~ s/\&[Ll]ambda;/lambda/g;
    $bit =~ s/\&[Xx]i;/xi/g;

    $bit =~ s/\&plusmn;/\(+\/-\)/g;
    $bit =~ s/\&rarr;/->/g;
    $bit =~ s/\&prime;?/\'/g;
    $bit =~ s/\&deg;//g;
    $bit =~ s/\&acirc;?/y/g;
    $bit =~ s/\&mdash;/\-/g;
    $bit =~ s/\&\#8596;/\<\-\>/g;
    $bit =~ s/\&\#8217;?/\'/g;
    $bit =~ s/\&\#8242;/\'/g;
    $bit =~ s/\&\#8222;/\'/g;

    $bit =~ s/\<\/?[ib]\/?\>//ig;
    $bit =~ s/\<\/?em\>//g;
    $bit =~ s/\<\/?small\>//g;
    $bit =~ s/\<\/?su[pb]\>//ig;
    $bit =~ s/\<\/?span\s?.*?\>//g;
    $bit =~ s/\<a\shref.*?\>//g;
    $bit =~ s/\<\/a\>//g;

    if($bit =~ /;/){
	print TRERR "WARNING: semi-colon found in name: ",$bit,"\n";
	$bit =~ s/;/,/g;
    }
    
    if($bit =~ /\&/){
	print TRERR "WARNING: ampersand found in name: ",$bit,"\n";
    }

    return $bit;
}

my $Database_Root = "/homes/seaver/Biochemistry_Mirrors/bioinformatics.ai.sri.com/metacyc/dist/22.5/data/";
my $Structure_Root = "/homes/seaver/Projects/ModelSEEDDatabase/Biochemistry/Structures/MetaCyc/";
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

my %CycSEED = ();
open(FH, "< ".$ENV{SEAVER_PROJECT}."Cyc/Meta/Meta_SEED_Compartments_v3.txt");
while(<FH>){
    chomp;
    @temp=split /\t/;
#    $CycSEED{$temp[0]}=$temp[1];
    $CycSEED{$temp[0]}=$temp[0];
}
close(FH);

#if($Dir ne "Meta" && -f $Dir."/".$Dir."_SEED_Compartments.txt"){
#    open(FH, "< ".$Dir."/".$Dir."_SEED_Compartments.txt");
#    while(<FH>){
#	chomp;
#	@temp=split /\t/;
#	$CycSEED{$temp[0]}=$temp[1];
#    }
#    close(FH);
#}

my $Compounds={};
open(FH, "< ".$Database_Root."/compounds.dat");
open(KOUT, "> MetaCyc_Output/KEGG_Compound_Match.txt");
open(TRERR, "> MetaCyc_Output/TranslateHTML_Warnings.txt");
open(IDEN, "> MetaCyc_Output/IdenticalNames_Warnings.txt");
my ($ID,$Formula,$Mass,$Charge)=("","",10000000,10000000);
my ($Field,$Data)=("","");
my %Names=();
my %Names_ID=();
my %DBLinks=();
my ($DB,$Link)=("","");
my %Cpd_Rxns=();
my %Radicals=();
while(<FH>){
    chomp;
    if(substr($_,0,1) eq "#"){next;}
    @temp=split/\s-\s/;

    $Field=$temp[0];
    $Data=join("-",@temp[1..$#temp]);
    
    $Field=~s/\s*$//;
    $Data=~s/^\s*//;
    $Data=~s/\s*$//;

    if($Field eq "UNIQUE-ID"){
	#07/25
	#All UNIQUE-IDs are UNIQUE
	$Data =~ s/\s/_/g;
	$ID=$Data;
    }elsif($Field eq "CHEMICAL-FORMULA"){
	#remove parenthesis
	$Data =~ s/[\(\)]//g;
	my ($atom,$stoich) = split(/\s/,$Data);
	$atom = substr($atom,0,1).lc(substr($atom,1,1));

	#remove 1 in cases of singles
	#Reminder CHEMICAL-FORMULA does different
	#atoms on different lines
	$Data = $atom;
	$Data .= $stoich > 1 ? $stoich : "";
	$Formula.=$Data;

    }elsif($Field eq "SYNONYMS" || $Field =~ /NAME$/){
	#07/25
	#metacyc counts
	#1 N-1-NAME
	#6 N-NAME
	#27 N+1-NAME
	#46 ABBREV-NAME
	#138 SYSTEMATIC-NAME
	#9187 COMMON-NAME
	next if $Field eq "N-NAME" || $Field =~ /^N[+-]1-NAME$/;
	next if !$Data;

	$Data=translate_html_entities($Data);
	if(exists($Names_ID{$Data}) && !exists($Names_ID{$Data}{$ID})){
            print IDEN "Skipping identical name $Data for $ID; already found in ".join("|",keys %{$Names_ID{$Data}})."\n";
	}else{
	    $Names{$Data}=1;
	}
	$Names_ID{$Data}{$ID}=1;
    }elsif($Field eq "CHARGE"){
	$Charge=$Data;
    }elsif($Field eq "MOLECULAR-WEIGHT"){
	#rare cases of trailing period
	$Data =~ s/\.$//;
	$Mass=$Data;
    }elsif($Field eq "DBLINKS"){
	#remove parenthesis
	$Data =~ s/[\(\)]//g;

	@temp=split /\s/,$Data;
	($DB,$Link)=($temp[0],$temp[1]);

	#remove bars
	$DB =~ s/\|//g;
	#remove \"
	$Link =~ s/\"//g;

	if($DB eq "Wikipedia"){
	    $Names{$Link}=1;
	}
	#remove "-CPD"
	$DB =~ s/-CPD$//;
	if($DB eq "LIGAND"){
	    $Link =~ s/^c/C/;
	    print KOUT $ID,"\t",$Link,"\n";
	#    $KEGG=$Link;
	}
	$DBLinks{$DB.":".$Link}=1;
    }elsif($Field eq "TYPES"){
	if($Data eq "Radicals"){
	    $Radicals{$ID}=1;
	}
    }elsif(substr($_,0,2) eq "//"){

	$Compounds->{$ID}={METACYC=>$ID,NAMES=>"",FORMULA=>"null",CHARGE=>10000000,MASS=>10000000};
	$Cpd_Rxns{$ID}=();

	$Compounds->{$ID}{NAMES}=join("|",sort keys %Names);
	$Compounds->{$ID}{NAMES}=$ID if !$Compounds->{$ID}{NAMES};
	$Compounds->{$ID}{FORMULA}=$Formula if $Formula && $Formula ne "null";
	$Compounds->{$ID}{CHARGE}=$Charge;
	$Compounds->{$ID}{MASS}=$Mass;

	#reset records
	($ID,$Formula,$Mass,$Charge)=("","",10000000,10000000);
	undef(%Names);
	undef(%DBLinks);
    }
}
close(FH);
close(KOUT);
close(IDEN);

open(RAD, "> MetaCyc_Output/Radicals.txt");
print RAD join("\n",sort keys %Radicals),"\n";
close(RAD);

#Find classes
open(FH, "< ".$Database_Root."/classes.dat");
open(CLASS, "> MetaCyc_Output/Classes_not_Compounds.dat");
$ID="";
my $Classes={};
%Names=();
($Field,$Data)=("","");
while(<FH>){
    chomp;
    if(substr($_,0,1) eq "#"){next;}
    @temp=split/\s-\s/;

    $Field=$temp[0];
    $Data=join("-",@temp[1..$#temp]);
    
    $Field=~s/\s*$//;
    $Data=~s/^\s*//;
    $Data=~s/\s*$//;

    if($Field eq "UNIQUE-ID"){
	$Data =~ s/\s/_/g;
	$ID=$Data;
    }elsif($Field eq "TYPES"){
	#$Classes->{$ID}=$Data;
    }elsif($Field eq "SYNONYMS" || $Field =~ /NAME$/){
	next if $Field eq "N-NAME" || $Field =~ /^N[+-]1-NAME$/;
	next if !$Data;

	$Data=translate_html_entities($Data);
	$Names{$Data}=1;
    }elsif(substr($_,0,2) eq "//"){

	$Classes->{$ID}={METACYC=>$ID,NAMES=>$ID,FORMULA=>"null",CHARGE=>10000000,"MASS"=>10000000};
	$Compounds->{$ID}{NAMES}=$ID if exists($Compounds->{$ID});
	print CLASS $ID,"\n" if !exists($Compounds->{$ID});

	#special condition accounting for acyl-carrier-proteins
	if($ID =~ /ACPs?/i){
	    $Names{$ID}=1;
	    my $Name_String=join("|",sort keys %Names);
	    $Classes->{$ID}{NAMES}=$Name_String;
	    $Compounds->{$ID}{NAMES}=$Name_String if exists($Compounds->{$ID});
	}

	$Cpd_Rxns{$ID}=();

	$ID="";
	undef(%Names);
    }    
}

#Add electrons
$Classes->{"E-"}={METACYC=>"E-",NAMES=>"E-",FORMULA=>"null",UNCHARGED_FORMULA=>"null",CHARGE=>-1,"MASS"=>0,STRINGCODE=>"nostringcode",GROUPS=>"nogroups"};
close(CLASS);

#Get gene references from genes.dat
open(FH, "< ".$Database_Root."/genes.dat");
($Field,$Data)=("","");
my ($peg,$accession)=("","");
my @product=();
my %Genes=();

while(<FH>){
    chomp;
    if(substr($_,0,1) eq "#"){next;}
    if(substr($_,0,1) eq "/" && substr($_,0,2) ne "//"){next;}

    @temp=split/\s-\s/;

    $Field=$temp[0];
    $Data=join("-",@temp[1..$#temp]);
    $Field=~s/\s*$//;
    $Data=~s/^\s*//;
    $Data=~s/\s*$//;

    if($Field eq "UNIQUE-ID"){
	$Data =~ s/\s/_/g;
	$ID=$Data;
    }elsif($Field eq "PRODUCT"){
	push(@product,$Data);
    }elsif($Field eq "ACCESSION-1"){
	$accession=uc($Data);
    }elsif($Field eq "//"){
	if($accession eq ""){$accession = $ID;}
	foreach my $product (@product){
	    if(exists($Genes{$product})){
		if((split/-/,$product)[0] eq $ID){
		    $Genes{$product}{"GENE"}{$accession}=1;
		    $Genes{$product}{"ID"}=$ID;
		}
	    }else{
		$Genes{$product}{"GENE"}{$accession}=1;
		$Genes{$product}{"ID"}=$ID;
	    }
	}
	($ID,$accession)=("","");
	undef(@product);
    }
}
close(FH);

#Get compartments from proteins.dat
open(FH, "< ".$Database_Root."/proteins.dat");
($Field,$Data)=("","");
my $gene="";
my @components=();
my $Proteins={};
my %Prots=();

my %GlobalCmpts=();

open(WARN, "> MetaCyc_Output/Gene_Protein_Inconsistencies.txt");
while(<FH>){
    chomp;
    if(substr($_,0,1) eq "#"){next;}
    if(substr($_,0,1) eq "/" && substr($_,0,2) ne "//"){next;}

    @temp=split/\s-\s/;

    $Field=$temp[0];
    $Data=join("-",@temp[1..$#temp]);
    $Field=~s/\s*$//;
    $Data=~s/^\s*//;
    $Data=~s/\s*$//;

    if($Field eq "UNIQUE-ID"){
	$Data =~ s/\s/_/g;
	$ID=$Data;
    }elsif($Field eq "LOCATIONS"){
	$Prots{$ID}{"Loc"}{$Data}=1 unless $Data eq "NIL";
	$GlobalCmpts{$Data}=1;
    }elsif($Field eq "GENE"){
	$gene=$Data;
    }elsif($Field eq "COMPONENTS"){
	push(@components,$Data);
    }elsif($Field eq "//"){
	$Proteins->{$ID}={METACYC=>$ID,NAMES=>$ID,FORMULA=>"null",CHARGE=>10000000,"MASS"=>10000000};

	if(exists($Genes{$ID}) && $Genes{$ID}{"ID"} ne $gene){
	    print WARN "WARNING (proteins.dat): Inconsistency in genes used: ",$ID,"\t",$gene,"\t",$Genes{$ID}{"ID"},"\n";
	}

	if(!exists($Genes{$ID})){
	    foreach my $c(@components){
		if(exists($Genes{$c})){
		    foreach my $gene (sort keys %{$Genes{$c}{"GENE"}}){
			$Prots{$ID}{"Genes"}{$gene}=1;
		    }
		}else{
		    print WARN "WARNING (proteins.dat): Protein component not found in Genes.dat :",$c,"\n";
		}
	    }
	}else{
	    foreach my $gene (sort keys %{$Genes{$ID}{"GENE"}}){
		$Prots{$ID}{"Genes"}{$gene}=1;
	    }
	}
	$gene="";
	undef(@components);
	$ID="";
    }
}
close(FH);
close(WARN);

#Get reaction direction from enzrxns.dat
#Also link proteins to reactions
open(FH, "< ".$Database_Root."/enzrxns.dat");
open(WARN, "> MetaCyc_Output/Enzyme_Protein_Inconsitencies.txt");
($Field,$Data)=("","");
my ($Rev,$Enz,$Rxn)=("","","");
my %Enzs=();
my %Rxns=();
my %EnzRxns=();
while(<FH>){
    chomp;
    if(substr($_,0,1) eq "#"){next;}
    if(substr($_,0,1) eq "/" && substr($_,0,2) ne "//"){next;}

    @temp=split/\s-\s/;

    $Field=$temp[0];
    $Data=join("-",@temp[1..$#temp]);
    $Field=~s/\s*$//;
    $Data=~s/^\s*//;
    $Data=~s/\s*$//;

    if($Field eq "UNIQUE-ID"){
	$Data =~ s/\s/_/g;
	$ID=$Data;
    }elsif($Field eq "ENZYME"){
	$Enz=$Data;
    }elsif($Field eq "REACTION"){
	$Rxn=$Data;
    }elsif($Field eq "REACTION-DIRECTION"){
	#08/26/11
	#Only ONE REACTION-DIRECTION attribute per record
	$Rev=$Data;
    }elsif(substr($_,0,2) eq "//"){
	if($Enz ne "" && !exists($Prots{$Enz})){
	    print WARN "WARNING (enzrxns.dat): Enzyme not found in proteins.dat :",$ID,"\t",$Enz,"\n";
	}

	$Enzs{$Enz}{$Rev}=1 if $Enz && $Rev;
	$Rxns{$Rxn}{$Rev}=1 if $Rxn && $Rev;
	$EnzRxns{$ID}={"Enz"=>$Enz,"Rxn"=>$Rxn};

	($ID,$Enz,$Rxn,$Rev)=("","","","");
    }
}
close(FH);
close(WARN);

open(FH, "< ".$Database_Root."/reactions.dat");
open(INCHIWARN, "> MetaCyc_Output/InChI_Warnings.txt");
open(RXNWARN, "> MetaCyc_Output/Reaction_Inconsistencies.txt");
open(RXN, "> MetaCyc_Output/Reactions.txt");
open(TRAN, "> MetaCyc_Output/Transporters.txt");
open(KOUT, "> MetaCyc_Output/KEGG_Reaction_Match.txt");
open(GPEGS, "> MetaCyc_Output/EC_Genes_Pegs.txt");
$Rev="";
my $Rxn_ID="";
my @EC=();
my @Left=();
my @Right=();
%Names=();
my $Loc="";
($Field,$Data)=("","");
my ($Eq_Side,$Eq_Count)=("",0);
my $In_Pathway=0;
my %Reactions=();
my @ID_Count=();
my %COMP_COUNT=();
my %EnzRxn_Links=();
my %RevList=();
my %CmpList=();
my %TrnsList=();
my %MissCmpts=();
my %PegList=();
my %TypesList=();
my %OrphanCs=();
my %UnCompounds=();
my $SpontRxn="N";
my %CompRxns=();
undef(%DBLinks);

while(<FH>){
    chomp;
    if(substr($_,0,1) eq "#"){next;}
    if(substr($_,0,1) eq "/" && substr($_,0,2) ne "//"){next;}

    @temp=split/\s-\s/;

    $Field=$temp[0];
    $Data=join("-",@temp[1..$#temp]);
    $Field=~s/\s*$//;
    $Data=~s/^\s*//;
    $Data=~s/\s*$//;

    #remove bars
    $Data =~ s/^\|//;
    $Data =~ s/\|$//;

    if($Field eq "UNIQUE-ID"){
	#07/26
	#All UNIQUE-IDs are UNIQUE
	$Data =~ s/\s/_/g;
	$ID=$Data;
	push(@ID_Count,$ID) if exists($Compounds->{$ID});
    }elsif($Field eq "SYNONYMS" || $Field =~ /NAME$/){
	#07/26
	#metacyc counts
	#3913 SYSTEMATIC-NAME
	#5874 COMMON-NAME
	#17727 SYNONYMS

	$Data=translate_html_entities($Data);
	$Names{$Data}=1;
    }elsif($Field eq "TYPES"){
	next if $Data =~ /^EC-\d+/;
	$TypesList{$Data}=1;
	if($Data eq "Composite-Reactions"){
	    $CompRxns{$ID}=1;
	}
    }elsif($Field eq "REACTION-LIST"){
	#For removing half-reactions
	#but balancing code should make these redundant	
	#$RXNList{$ID}{$Data}=1;
    }elsif($Field eq "LEFT"){
	$Data =~ s/\s/_/g;
	push(@Left,{"ID"=>$Data,"COEFF"=>1,"COMP"=>""});
	$Cpd_Rxns{$Data}{$ID}=1;

	$Eq_Side=$Field;
	$Eq_Count=$#Left;

	if(!exists($Compounds->{$Data})){
	    $OrphanCs{$Data}{$ID}=1;
	}
    }elsif($Field eq "RIGHT"){
	$Data =~ s/\s/_/g;
	push(@Right,{"ID"=>$Data,"COEFF"=>1,"COMP"=>""});
	$Cpd_Rxns{$Data}{$ID}=1;

	$Eq_Side=$Field;
	$Eq_Count=$#Right;
	
	if(!exists($Compounds->{$Data})){
	    $OrphanCs{$Data}{$ID}=1;
	}
    }elsif($Field eq "^COMPARTMENT"){
	next unless substr($Data,0,3) eq "CCO";
	$GlobalCmpts{$Data}=1;
	if($Eq_Side eq "LEFT"){
	    $Left[$Eq_Count]->{"COMP"}=$Data;
	}elsif($Eq_Side eq "RIGHT"){
	    $Right[$Eq_Count]->{"COMP"}=$Data;
	}
	$COMP_COUNT{$Data}=1;
    }elsif($Field eq "^COEFFICIENT"){
	$Data=translate_html_entities($Data);
	$Data =~ s/(\d)n/$1/;
	if($Data =~ /^\D{1}$/){$Data="1";}
	if($Data =~ /\D+\D/){$Data="2";}

	if($Eq_Side eq "LEFT"){
	    $Left[$Eq_Count]->{"COEFF"}=$Data;
	}elsif($Eq_Side eq "RIGHT"){
	    $Right[$Eq_Count]->{"COEFF"}=$Data;
	}
    }elsif($Field eq "ENZYMATIC-REACTION"){
	$EnzRxn_Links{$Data}=1;
    }elsif($Field eq "EC-NUMBER"){
	push(@EC,$Data);
    }elsif($Field eq "REACTION-DIRECTION"){
	$Rev=$Data;
    }elsif($Field eq "RXN-LOCATIONS"){
	if($Data ne "NIL"){
	    $Loc=$Data;
	}
    }elsif($Field eq "SPONTANEOUS?"){
	$SpontRxn="Y" if $Data eq "T";
    }elsif($Field eq "IN-PATHWAY"){
	$In_Pathway=1;
    }elsif($Field eq "DBLINKS"){
	#remove parenthesis
	$Data =~ s/[\(\)]//g;

	@temp=split /\s/,$Data;
	($DB,$Link)=($temp[0],$temp[1]);

	#remove bars
	$DB =~ s/\|//g;
	#remove \"
	$Link =~ s/\"//g;

	#remove "-RXN"
	$DB =~ s/-RXN$//;
	if($DB eq "LIGAND"){
	    print KOUT $ID,"\t",$Link,"\n";
	#    $KEGG=$Link;
	}
	$DBLinks{$DB.":".$Link}=1;
    }elsif(substr($_,0,2) eq "//"){
	if(scalar(@Left)==0 || scalar(@Right)==0){
	    print RXNWARN "WARNING (reactions.dat) no reactants and/or products: ",$ID,"\t",scalar(@Left),"\t",scalar(@Right),"\n";
	    ($ID,$Rev,$Loc,$In_Pathway,$SpontRxn)=("","","",0,"N");
	    undef(@EC);
	    undef(%Names);
	    undef(@Left);
	    undef(@Right);
	    undef(%COMP_COUNT);
	    undef(%EnzRxn_Links);
	    undef(%TypesList);
	    undef(%DBLinks);
	    next;
	}

	undef(%PegList);
	my %genes=();
	foreach my $k(keys %EnzRxn_Links){
	    foreach my $gene (keys %{$Prots{$EnzRxns{$k}{"Enz"}}{"Genes"}}){
		$PegList{$gene}=1;
		$genes{$gene}=1 if $gene !~ /peg/;
	    }
	}

	if($SpontRxn eq "Y"){
	    if(scalar(keys %PegList)>0){
		print RXNWARN "WARNING (reactions.dat) spontaneous reaction $ID has associated pegs\n";
	    }
	    $PegList{"Spont"}=1;
	}

	print GPEGS $ID,"\t",join("|", @EC),"\t",join("|",keys %genes),"\n";

	undef(%CmpList);
	undef(%RevList);
	foreach my $k(keys %EnzRxn_Links){
	    if(exists($Enzs{$EnzRxns{$k}{"Enz"}})){
		foreach my $k2 (keys %{$Enzs{$EnzRxns{$k}{"Enz"}}}){
		    $RevList{$k2}=1 if $k2 ne "REVERSIBLE";
		}
	    }
	    if(exists($Rxns{$EnzRxns{$k}{"Rxn"}})){
		foreach my $k2 (keys %{$Rxns{$EnzRxns{$k}{"Rxn"}}}){
		    $RevList{$k2}=1 if $k2 ne "REVERSIBLE"; 
		}
	    }
	    foreach my $k2(grep { $_ ne "" } keys %{$Prots{$EnzRxns{$k}{"Enz"}}{"Loc"}}){
		if((() = $k2 =~ /CCO/g)==2){
		    $k2 =~ /^(CCO.*?)-(CCO.*?)$/;
		    my @locs=($1,$2);
		    foreach my $loc (@locs){
			$CmpList{$loc}=1 if exists($CycSEED{$loc});
			$MissCmpts{$loc}=1 if !exists($CycSEED{$loc});
		    }
		}else{
		    $CmpList{$k2}=1 if exists($CycSEED{$k2});
		    $MissCmpts{$k2}=1 if !exists($CycSEED{$k2});
		}
	    }
	}

	if($Loc){
	    if((() = $Loc =~ /CCO/g)==2){
		$Loc =~ /^(CCO.*?)-(CCO.*?)$/;
		my @locs=($1,$2);
		foreach my $loc (@locs){		    	
		    $GlobalCmpts{$loc}=1;
		    $CmpList{$loc}=1 if exists($CycSEED{$loc});
		    $MissCmpts{$loc}=1 if !exists($CycSEED{$loc}) && $Loc ne "";
		}
	    }else{
		$GlobalCmpts{$Loc}=1;
		$CmpList{$Loc}=1 if exists($CycSEED{$Loc});
		$MissCmpts{$Loc}=1 if !exists($CycSEED{$Loc}) && $Loc ne "";
	    }
	}elsif(scalar(keys %CmpList)==0){
	    $CmpList{"CCO-CYTOSOL"}=1;
	}

	undef(%TrnsList);
	if((() = $Loc =~ /CCO/g)==2){
	    $Loc =~ /^(CCO.*?)-(CCO.*?)$/;
	    my @locs=($1,$2);
	    my $yes_transport=1;
	    my %tmp=();
	    foreach my $l (@locs){
		if(!exists($CycSEED{$l}) && $l ne "CCO-OUT" && $l ne "CCO-IN"){
		    $yes_transport=0;
		    last;
		}else{
		    if(exists($CycSEED{$l})){
			$tmp{$CycSEED{$l}}=1;
		    }else{
			$tmp{$l}=1;
		    }
		}
	    }

	    $yes_transport=0 if scalar(keys %tmp) < 2;

	    if($yes_transport){
		my %trans=();
		my $innermost="";
		for(my $i=0;$i<2;$i++){
		    $locs[$i] = "CCO-EXTRACELLULAR" if $locs[$i] eq "CCO-OUT";
		    $locs[$i] = "CCO-CYTOSOL" if $locs[$i] eq "CCO-IN" || $locs[$i] eq "CCO-UNKNOWN-SPACE";
		    if(!exists($trans{"IN"}) && ($locs[$i] eq "CCO-CYTOSOL" || $locs[$i] eq "CCO-MIT-LUM" || $locs[$i] eq "CCO-CHL-THY-LUM")){
			$trans{"IN"}=$locs[$i];
		    }
		    if($locs[$i] eq "CCO-EXTRACELLULAR"){
			$trans{"OUT"}=$locs[$i];
		    }
		    
		    if($locs[$i] ne "CCO-EXTRACELLULAR" && $locs[$i] ne "CCO-CYTOSOL"){
			$innermost=$locs[$i];
		    }
		}

		#Unusual case of two compartments being same
		if($locs[0] eq $locs[1]){
		    $innermost="";
		    undef(%trans);
		    %CmpList=($locs[0]=>1);
		}else{

		    for(my $i=0;$i<2;$i++){
			next if $locs[$i] eq "CCO-EXTRACELLULAR" || $locs[$i] eq "CCO-CYTOSOL";
		    
			if(exists($trans{"IN"})){
			    $trans{"OUT"}=$locs[$i];
			}elsif(exists($trans{"OUT"})){
			    $trans{"IN"}=$locs[$i];
			}
		    }

		    $innermost = "CCO-CYTOSOL" if $innermost eq "";
		    $TrnsList{$innermost}{"IN"} = $trans{"IN"};
		    $TrnsList{$innermost}{"OUT"} = $trans{"OUT"};
		    undef(%CmpList);
		}
	    }else{
		#As of November 2014, this has happened once, because two separate Cyc compartments had the same SEED compartment
		print TRAN "Locs (",join("|",@locs),") translate to same SEED Cmpts: ",join("|",keys %tmp),".\n";
		if(scalar(keys %tmp)==0){
		    %CmpList=("CCO-CYTOSOL"=>1);
		}else{
		    %CmpList=($locs[0]=>1);
		}
	    }
	}elsif(scalar(keys %COMP_COUNT)>1){

	    if(scalar(keys %CmpList)==1 && (keys %CmpList)[0] eq "CCO-CYTOSOL"){
		    $TrnsList{"CCO-CYTOSOL"}{"IN"}="CCO-CYTOSOL";
		    $TrnsList{"CCO-CYTOSOL"}{"OUT"}="CCO-EXTRACELLULAR";
	    }else{
		#only transporter if SEED compartment (all Cyc databases bar Plant/Meta should have all)
		foreach my $c ( grep { exists($CycSEED{$_}) && $CycSEED{$_} ne 'c' } keys %CmpList){
		    my ($OUT,$IN) = ("CCO-CYTOSOL",$c);
		    if($CycSEED{$c} eq 'e' || $CycSEED{$c} eq 'p' || $CycSEED{$c} eq 'w'){
			($OUT,$IN)=($c,"CCO-CYTOSOL");
		    }
		    $TrnsList{$c}{"IN"}=$IN;
		    $TrnsList{$c}{"OUT"}=$OUT;
		    
		    if($IN eq $OUT){
			print TRAN "Warning: Transporter reaction $ID happening within one compartment: $IN\n";
		    }
		}
	    }
	}

	#update RevList so that it has LEFT-RIGHT or RIGHT-LEFT only
	my %tmpList=%RevList;
	undef(%RevList);
	foreach my $r (keys %tmpList){
	    $RevList{substr($r,-13)}=1;
	}
	$Rev=substr($Rev,-13);	    

	if(scalar(keys(%RevList))>1){
	    print RXNWARN "WARNING (enzrxns.dat) Multiple Directions listed for ",$ID," : ",join("|",keys(%RevList)),"\n";
	    $Rev="REVERSIBLE";
	    #print "Multiple Revs\t",$Rev,"\n";
	}elsif(scalar(keys(%RevList))==1 && (keys %RevList)[0] ne "REVERSIBLE"){
	    if((keys(%RevList))[0] ne $Rev && $Rev ne "REVERSIBLE" && $Rev ne ""){
		print RXNWARN "WARNING (enzrxns.dat/reactions.dat) Multiple Directions listed for ",$ID," : ",join("|",((keys(%RevList))[0],$Rev)),"\n";
	    }
	    $Rev=(keys %RevList)[0];
	    #print "ENZRXN Revs\t",$Rev,"\n";
	}elsif(!$In_Pathway && $Rev ne ""){
	    #print "REACTION Revs\t",$Rev,"\t",$In_Pathway,"\n";
	}else{
	    $Rev="REVERSIBLE";
	    #print "Default Revs\t",$Rev,"\n";
	}

	if(scalar(keys %TrnsList)){
	    my %tmpList=();
	    
	    foreach my $c (keys %TrnsList){
		if(!defined($TrnsList{$c}{"IN"})){
		    print $ID,"\t",$c,"\n";
		}
		$tmpList{$TrnsList{$c}{"IN"}}=1;
		$tmpList{$TrnsList{$c}{"OUT"}}=1;
		$TrnsList{$c}{"IN"}=$CycSEED{$TrnsList{$c}{"IN"}};
		$TrnsList{$c}{"OUT"}=$CycSEED{$TrnsList{$c}{"OUT"}};

		generate_reaction(\@Left, \@Right,[$TrnsList{$c}{"IN"},$TrnsList{$c}{"OUT"}],$Rev,1);
		print TRAN "Transporter found (reactions.dat): ",$ID,"\t",join("|",keys %tmpList),"\t",join("|",($TrnsList{$c}{"IN"},$TrnsList{$c}{"OUT"})),"\n";
	    }
	}else{
	    if(!scalar(keys %CmpList)){
		$CmpList{"CCO-CYTOSOL"}=1;
	    }

	    my %tmpList=%CmpList;
	    undef(%CmpList);
	    foreach my $c (keys %tmpList){
		$CmpList{$CycSEED{$c}}=1;
	    }

	    print RXN "Reaction found (reactions.dat): ",$ID,"\t",join("|",keys %tmpList),"\t",join("|",keys %CmpList),"\n";

	    foreach my $c (keys %CmpList){
		generate_reaction(\@Left,\@Right,[$c],$Rev,1);
	    }
	}

	($ID,$Rev,$Loc,$In_Pathway,$SpontRxn)=("","","",0,"N");
	undef(@EC);
	undef(%Names);
	undef(@Left);
	undef(@Right);
	undef(%COMP_COUNT);
	undef(%EnzRxn_Links);
	undef(%TypesList);
	undef(%DBLinks);
    }
}
close(FH);
close(INCHIWARN);
close(RXNWARN);
close(TRAN);
close(RXN);
close(KOUT);
close(GPEGS);

open(LOC, "> MetaCyc_Output/Protein_Locations.txt");
foreach my $ID ( grep { exists($Prots{$_}) && exists($Prots{$_}{'Loc'}) && scalar(keys %{$Prots{$_}{'Loc'}})>0 } keys %Prots ){
    print STDERR "Warning: Missing compartments: ".join("|", grep { !exists($CycSEED{$_}) } keys %{$Prots{$ID}{'Loc'}}),"\n" if scalar( grep { !exists($CycSEED{$_}) } keys %{$Prots{$ID}{'Loc'}})>0;
    print LOC $ID,"\t",(exists($Prots{$ID}{"Genes"}) ? join("|", grep { $_ !~ /peg\./ } keys %{$Prots{$ID}{"Genes"}}) : ""),"\t",join("|", sort keys %{$Prots{$ID}{'Loc'}}),"\t",join("|", sort map { $CycSEED{$_} } keys %{$Prots{$ID}{'Loc'}}),"\n";
}
close(LOC);

open(COMP, "> MetaCyc_Output/Compartments.txt");
foreach my $c (keys %GlobalCmpts){
    print COMP $c,"\n";
}
close(COMP);

open(COMP, "> MetaCyc_Output/Missing_Compartments.txt");
foreach my $c (keys %MissCmpts){
    print COMP $c,"\n";
}
close(COMP);

open(CMPRXN, "> MetaCyc_Output/Composite_Reactions.txt");
foreach my $r (keys %CompRxns){
    print CMPRXN $r,"\n";
}
close(CMPRXN);

my %Remove_Names=();
#my %Remove_Names=("phosphoglycerate"=>1,"aldehyde"=>1);
open(OUT, "> ".$Output_Root."MetaCyc_Compounds.tbl");
my @Headers=("ID","NAMES","FORMULA","CHARGE","MASS");
print OUT join("\t",@Headers)."\tInChI\tSMILE\n";
foreach my $cpd (keys %$Compounds){
    print OUT $cpd,"\t";
    foreach my $h ( grep { $_ ne "ID" } @Headers){
	if($h eq "NAMES"){
	    $Compounds->{$cpd}{$h} = join("|", grep { !exists($Remove_Names{$_}) && length($_)>1 } split(/\|/,$Compounds->{$cpd}{$h}));
	}

	print OUT $Compounds->{$cpd}{$h};
	print OUT "\t"; # unless $h eq $Headers[$#Headers];
    }
    foreach my $type ("InChI","SMILE"){
	if(exists($Structures{$cpd}) && exists($Structures{$cpd}{$type}) && $Structures{$cpd}{$type} ne ""){
	    print OUT $Structures{$cpd}{$type};
	}else{
	    print OUT "-";
	}
	print OUT "\t";
    }
    print OUT "\n";
}

open(MAP, "> MetaCyc_Output/Compound_Reaction_Mapping.txt");
foreach my $cpd (sort keys %Cpd_Rxns){
    print MAP $cpd,"\t",join("|",sort keys %{$Cpd_Rxns{$cpd}}),"\n";
}
close(MAP);

my %Missing_Rxns=();
open(MSRXN, "> MetaCyc_Output/Excluded_Reactions_Missing_Compounds.txt");
open(ORPH, "> MetaCyc_Output/NonCompound_Reactions.txt");
foreach my $cpd (keys %OrphanCs){
    if(exists($Classes->{$cpd})){
	print ORPH "Class\t",$cpd,"\t",join("|",sort keys %{$OrphanCs{$cpd}}),"\n";
    
	print OUT $cpd,"\t";
	foreach my $h ( grep { $_ ne "ID" } @Headers){
	    if(!exists($Classes->{$cpd}{$h})){
		print "Missing $h\n";
	    }
	    print OUT $Classes->{$cpd}{$h};
	    print OUT "\t";
	}
	foreach my $type ("InChI","SMILE"){
	    if(exists($Structures{$cpd}) && exists($Structures{$cpd}{$type}) && $Structures{$cpd}{$type} ne ""){
		print OUT $Structures{$cpd}{$type};
	    }else{
		print OUT "-";
	    }
	    print OUT "\t";
	}
	print OUT "\n";
    }elsif(exists($Proteins->{$cpd})){
	print ORPH "Protein\t",$cpd,"\t",join("|",sort keys %{$OrphanCs{$cpd}}),"\n";
	print OUT $cpd,"\t";
	foreach my $h ( grep { $_ ne "ID" } @Headers){
	    print OUT $Proteins->{$cpd}{$h};
	    print OUT "\t";
	}
	foreach my $type ("InChI","SMILE"){
	    if(exists($Structures{$cpd}) && exists($Structures{$cpd}{$type}) && $Structures{$cpd}{$type} ne ""){
		print OUT $Structures{$cpd}{$type};
	    }else{
		print OUT "-";
	    }
	    print OUT "\t";
	}
	print OUT "\n";
    }else{
	foreach my $rxn(keys %{$OrphanCs{$cpd}}){
	    $Missing_Rxns{$rxn}=1;
	    print MSRXN "Missing: $cpd\t$rxn\n";
	}
    }
}
close(OUT);
close(ORPH);

open(OUT, "> ".$Output_Root."MetaCyc_Reactions.tbl");
print OUT "ID\tNAMES\tMETACYC\tEQUATION\tECs\n";
foreach my $rxn ( grep { !exists($Missing_Rxns{$_}) } keys %Reactions){
#    $Reactions{$rxn}{"Compartment"}="c" if !$Reactions{$rxn}{"Compartment"};
    print OUT $rxn,"\t",join("|",keys %{$Reactions{$rxn}{"Names"}}),"\t",$Reactions{$rxn}{"MetaCyc"},"\t",$Reactions{$rxn}{"Equation"},"\t",$Reactions{$rxn}{"EC"},"\n";
}
close(OUT);

sub generate_reaction{
    my $Left=shift;
    my $Right=shift;
    my $comp=shift;
    my $direction=shift;
    my $count=shift;
    my $tmpid=$ID;
#    $tmpid.=".".join('',sort @$comp) if $count;

    #print $tmpid,"\t",join("|", map { join("|", %$_) } @$Left),"\t",join("|", map { join("|", %$_) } @$Right),"--\t--",join("|",@$comp),"--\t--",$direction,"--\t--",$count,"--\n";

    $Reactions{$tmpid}={};
    $Reactions{$tmpid}{"Left"}=dclone($Left);
    assign_compartment($Reactions{$tmpid}{"Left"},@$comp);
    $Reactions{$tmpid}{"Right"}=dclone($Right);
    assign_compartment($Reactions{$tmpid}{"Right"},@$comp);

    foreach my $n(keys %Names){
	$Reactions{$tmpid}{"Names"}{$n}=1;
    }
    foreach my $db(keys %DBLinks){
	$Reactions{$tmpid}{"DBLinks"}{$db}=1;
    }
    
    $Reactions{$tmpid}{"EC"}=join("|",@EC);
    $Reactions{$tmpid}{"EC"}="-" if scalar(@EC)==0;
    ($Reactions{$tmpid}{"Equation"},$Reactions{$tmpid}{"Rev"},$Reactions{$tmpid}{"Compartment"})=cyc_reaction_2_string($Reactions{$tmpid}{"Left"},$Reactions{$tmpid}{"Right"},$direction);
#    print STDERR $ID,"\t",$tmpid,"\t",$comp,"\t",$Reactions{$tmpid}{"Equation"},"\n";

    $Reactions{$tmpid}{"PEGS"}=join("|",keys %PegList);
    $Reactions{$tmpid}{"MetaCyc"}=$ID;

#    print $tmpid,"\t",$Reactions{$tmpid}{"Equation"},"\n";
}

sub assign_compartment{
    my $array=shift;
    my @Comps=@_;

    if(scalar(@Comps)==1){
	foreach my $c (@$array){
#	    $c->{"COMP"}=$Comps[0];
	    $c->{"COMP"}="0";
	}
    }else{
	#MUST BE TWO!
	#$Comps[0] eq $IN
	#$Comps[1] eq $OUT
	foreach my $c(@$array){
	    if($c->{"COMP"} eq "CCO-OUT"){
#		$c->{"COMP"}=$Comps[1];
		$c->{"COMP"}="1";
	    }else{
#		$c->{"COMP"}=$Comps[0];
		$c->{"COMP"}="0";
	    }
	}
    }
}

sub cyc_reaction_2_string{
    my ($left,$right,$dir)=@_;

    my @ls=();
    my $string="";
    my $Strange="";

    my %Comps=();

    foreach my $l(@$left){
	$string="";
	if($l->{"COEFF"} =~ /\D/ && $l->{"COEFF"} !~ /\d\.\d/){
	    $Strange.=$l->{"COEFF"};
	}
	if($l->{"COEFF"}>1){
	    $string.="(".$l->{"COEFF"}.") ";
	}
	$string.=$l->{"ID"};
	$string.="[".$l->{"COMP"}."]" if $l->{"COMP"} ne "";
	$Comps{$l->{"COMP"}}=1;
	push @ls,$string;
    }

    my $eq=join(" + ",@ls);

    #Directionality
    my $rev="";
    if(substr($dir,-4) eq "LEFT"){
	$rev= " <= ";
    }elsif(substr($dir,-5) eq "RIGHT"){
	$rev=" => ";
    }else{
	$rev=" <=> ";
    }
    $eq.=$rev;

    my @rs=();
    $string="";
    foreach my $r(@$right){
	$string="";
	if($r->{"COEFF"} =~ /\D/ && $r->{"COEFF"} !~ /\d\.\d/){
	    $Strange.=$r->{"COEFF"};
	}
	if($r->{"COEFF"}>1){
	    $string.="(".$r->{"COEFF"}.") ";
	}
	$string.=$r->{"ID"};
	$string.="[".$r->{"COMP"}."]" if $r->{"COMP"} ne "";
	$Comps{$r->{"COMP"}}=1;
	push @rs,$string;
    }

    $eq.=join(" + ",@rs);

    my $reaction_comp= (scalar(keys %Comps)==1) ? (keys %Comps)[0] : "c";

    if($Strange){
	print STDERR "Unrecognized Coefficient ",$Strange," found in ",$ID,"\t",$eq,"\n";
    }

    return ($eq,$rev,$reaction_comp);
}

close(TRERR);
