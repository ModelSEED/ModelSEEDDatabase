package Pathways;
use warnings;
use strict;
my $header=1;
my @temp=();

my $MS_Path = "/vol/model-prod/";

my %Biochem_Sources = ("KEGG" => "ftp.bioinformatics.jp/kegg/pathway/",
		       "SRI" => "bioinformatics.ai.sri.com/",
		       "PMN" => "ftp.plantcyc.org/tmp/private/plantcyc/");

#Plans to include more BioCyc Databases, BiGG, and Subsystems
my %Biochem_Pathway_Files = ("KEGG" => $MS_Path.$Biochem_Sources{KEGG}."pathway",
			     "MetaCyc" => $MS_Path.$Biochem_Sources{SRI}."metacyc/pathways.dat",
			     "AraCyc" => $MS_Path.$Biochem_Sources{PMN}."aracyc/pathways.dat");

sub list_cyc_pathways{
    my @dirs=();
    foreach my $DB ("SRI","PMN"){
	opendir(my $dir, $MS_Path.$Biochem_Sources{$DB});
	push(@dirs, map { $DB.": ".$_ } grep { $_ =~ /cyc/ } readdir($dir));
	closedir($dir);
    }
    return @dirs;
}

sub get_pathways{
    my $hash = shift;

    if($hash->{db} eq "KEGG"){
	return get_kegg_pathways($hash);
    }elsif($hash->{db} eq "SRI" || $hash->{db} eq "PMN"){
	return get_cyc_pathways($hash);
    }else{
	return;
    }
}

sub get_kegg_pathways{
    my $hash = shift;

    my $Pathways;
    open(FH, "< ".$Biochem_Pathway_Files{$hash->{db}});
    my ($field,$data)=("","");
    my ($ENTRY,$CLASS,$NAME)=("","","");
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
	    $Pathways->{$ENTRY}{"Class"}=$CLASS;
	    $Pathways->{$ENTRY}{"Name"}=$NAME;
	    $Pathways->{$ENTRY}{"Reactions"}=[];
	    ($ENTRY,$CLASS,$NAME)=("","","");
	}
    }
    close(FH);

    my $Rxn_File = $MS_Path.$Biochem_Sources{KEGG}."../ligand/reaction/reaction";
    open(FH, "< $Rxn_File");
    undef(@temp);
    ($field,$data)=("","");
    $ENTRY="";
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
	if($field eq "PATHWAY"){
	    next if $data =~ /rn01100/; #skip generic "Metabolic pathways"
	    next if $data =~ /rn01110/; #skip generic "Biosynthesis of secondary metabolites"
	    next if $data =~ /rn01120/; #skip generic "Microbial metabolism in diverse environments"
	    @temp=split(/\s+/,$data);
	    $temp[0] =~ s/rn/map/;
	    push(@{$Pathways->{$temp[0]}{"Reactions"}},$ENTRY);
	}
    }
    close(FH);

    return $Pathways;
}


1;
