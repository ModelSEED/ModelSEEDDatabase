#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header=1;

######################
#Read in current roles
######################
open(FH, "< ../Mappings/probmodelseed-mapping/Mapping_Roles.txt");
my @Role_Headers = split(/\t/,<FH>);
chomp($Role_Headers[$#Role_Headers]);
my %Current_Role_Names = ();
my %Current_Role_IDs = ();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    $Current_Role_Names{$temp[1]}=$temp[0];
    $Current_Role_IDs{$temp[0]}=$temp[1];
}
close(FH);

######################
#Read in current complexes
######################
open(FH, "< ../Mappings/probmodelseed-mapping/Mapping_Complexes.txt");
my @Complex_Headers = split(/\t/,<FH>);
chomp($Complex_Headers[$#Complex_Headers]);
my %Current_Complexes = ();
while(<FH>){
    chomp;
     @temp=split(/\t/,$_,-1);
    $Current_Complexes{$temp[0]}{'name'}=$temp[1];
    $Current_Complexes{$temp[0]}{'roles'}=();
    foreach my $cpxrole (split(/\|/,$temp[2])){
	push(@{$Current_Complexes{$temp[0]}{'roles'}},$cpxrole);
    }
}
close(FH);

#"mscpx.54        cpx00054        role_ref:~/roles/id/msfr.1530;optionalRole:0;type:triggering;triggering:1"
#"complex_id      complex_name    complex_source  complex_type    role_id"

open(FH, "< ../Mappings/ComplexRoles.kegg.tsv");
my %Transform_IDs=();
my %New_Role_IDs = ();
my %New_Complexes=();
$header = 1;
while(<FH>){
    chomp;
    if($header){$header--;next}
    @temp=split(/\t/,$_,-1);

    my $role_id = $temp[4];

    my $Is_Present=0;
    if(exists($Current_Role_Names{$temp[5]})){
	$Transform_IDs{$temp[4]}=$Current_Role_Names{$temp[5]};
	$role_id = $Current_Role_Names{$temp[5]};
	$Is_Present=1;
    }

    #attempt aliases just in case
    foreach my $set_alias (split(/;/,$temp[8])){
	my ($set,$alias)=split(/:/,$set_alias);
	next if !$alias;

	if(exists($Current_Role_Names{$alias})){

	    #WARNING: Conflicting names and ids exist
	    if(exists($Transform_IDs{$temp[4]}) && $Transform_IDs{$temp[4]} ne $Current_Role_Names{$alias}){
		print "Warning: Conflicting ids/names\n";
		next;
	    }

	    $Transform_IDs{$temp[4]}=$Current_Role_Names{$alias};
	    $role_id = $Current_Role_Names{$alias};
	    $Is_Present=1;

	}
    }


    my $cpxrole = "role_ref:~/roles/id/".$role_id.";optionalRole:".$temp[12].";type:".$temp[10].";triggering:".$temp[11];

    if(!exists($New_Complexes{$temp[0]})){
	$New_Complexes{$temp[0]}{'name'}=$temp[1];
	$New_Complexes{$temp[0]}{'roles'}=[$cpxrole];
    }else{
	push(@{$New_Complexes{$temp[0]}{'roles'}},$cpxrole);
    }


    next if $Is_Present;
    $New_Role_IDs{$temp[4]}=$temp[5];
}
close(FH);

open(OUT, "> ../Mappings/probmodelseed-mapping/Mapping_Roles.txt");
print OUT join("\t",@Role_Headers),"\n";
foreach my $id (sort keys %Current_Role_IDs){
    print OUT $id,"\t",$Current_Role_IDs{$id},"\tnull\n";
}
foreach my $id (sort keys %New_Role_IDs){
    print OUT $id,"\t",$New_Role_IDs{$id},"\tnull\n";
}
close(OUT);

open(OUT, "> ../Mappings/probmodelseed-mapping/Mapping_Complexes.txt");
print OUT join("\t",@Complex_Headers),"\n";
foreach my $id (sort keys %Current_Complexes){
    print OUT $id,"\t",$Current_Complexes{$id}{'name'},"\t",join("|",sort @{$Current_Complexes{$id}{'roles'}}),"\n";
}
foreach my $id (sort keys %New_Complexes){
    print STDERR $id,"\n" if !exists($New_Complexes{$id}{'roles'});
    print OUT $id,"\t",$New_Complexes{$id}{'name'},"\t",join("|",sort @{$New_Complexes{$id}{'roles'}}),"\n";
}
close(OUT);
