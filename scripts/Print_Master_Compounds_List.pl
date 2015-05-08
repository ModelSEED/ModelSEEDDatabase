#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header=1;

my %Defaults = (charge => 0,
		isCofactor => 0,
		formula => "null");

#Load up the required modifications
open(FH, "< ../Biochemistry/compounds.master.mods");
my %Cpd_Mods=();
while(<FH>){
    chomp;
    @temp = split(/\t/,$_,-1);
    $Cpd_Mods{$temp[0]}{$temp[2]}=$temp[3];
}
close(FH);

#Start with original biochemistry
open(FH, "< ../Biochemistry/compounds.default.tsv");
my %Cpds=();
my @headers = split(/\t/,<FH>);
chomp $headers[$#headers];
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    for(my $i=1;$i<scalar(@temp);$i++){
	if(exists($Cpd_Mods{$temp[0]}) && exists($Cpd_Mods{$temp[0]}{$headers[$i]})){
	    $temp[$i] = $Cpd_Mods{$temp[0]}{$headers[$i]};
	}elsif(!$temp[$i] && exists($Defaults{$headers[$i]})){
	    $temp[$i] = $Defaults{$headers[$i]};
	}

	$Cpds{$temp[0]}{$headers[$i]}=$temp[$i];
    }
}
close(FH);

#Add new biochemistry
open(FH, "< ../Biochemistry/compounds.plantdefault.tsv");
$header=1;
while(<FH>){
    chomp;
    if($header){$header--;next;}
    @temp=split(/\t/,$_,-1);
    next if exists($Cpds{$temp[0]});

    for(my $i=1;$i<scalar(@temp);$i++){
	if(exists($Cpd_Mods{$temp[0]}) && exists($Cpd_Mods{$temp[0]}{$headers[$i]})){
	    $temp[$i] = $Cpd_Mods{$temp[0]}{$headers[$i]};
	}elsif(!$temp[$i] && exists($Defaults{$headers[$i]})){
	    $temp[$i] = $Defaults{$headers[$i]};
	}

	$Cpds{$temp[0]}{$headers[$i]}=$temp[$i];
    }
}
close(FH);

#Add aliases for KEGG and MetaCyc to match InChI structures
open(FH, "< ../Aliases/KEGG.aliases");
my %Aliases = ();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    foreach my $id (split(/\|/,$temp[2])){
	$Aliases{$id}{KEGG}=$temp[0];
    }
}
close(FH);

open(FH, "< ../Aliases/MetaCyc.aliases");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    foreach my $id (split(/\|/,$temp[2])){
	$Aliases{$id}{MetaCyc}=$temp[0];
    }
}
close(FH);

#Collect structures
open(FH, "< ../Structures/KEGG_Search_InChI.txt");
my %InChIs = ();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $InChIs{$temp[0]}=$temp[1];
}
close(FH);

open(FH, "< ../Structures/MetaCyc_Search_InChI.txt");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $InChIs{$temp[0]}=$temp[1];
}
close(FH);

#Print it all out
open(OUT, "> Master_Compound_List.tsv");
print OUT join("\t",@headers),"\n";
foreach my $cpd ( grep { $_ ne "cpd00000" } sort keys %Cpds){
    print OUT $cpd."\t";

    print OUT join("\t", map { $Cpds{$cpd}{$_} } grep { $_ ne "id" } @headers),"\t";

    #priortize InChIs from KEGG
    if(exists($Aliases{$cpd}{KEGG}) && exists($InChIs{$Aliases{$cpd}{KEGG}})){
	print OUT $InChIs{$Aliases{$cpd}{KEGG}}."\n";
    }elsif(exists($Aliases{$cpd}{MetaCyc}) && exists($InChIs{$Aliases{$cpd}{MetaCyc}})){
	print OUT $InChIs{$Aliases{$cpd}{MetaCyc}}."\n";
    }
}
