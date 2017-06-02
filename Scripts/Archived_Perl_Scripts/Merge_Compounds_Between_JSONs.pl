#!/usr/bin/env perl
use warnings;
use strict;
use JSON;

open(DB, "< ../Biochemistry/default_biochemistry.json");
my $Default="";
while(<DB>){
    $Default.=$_;
}
close(DB);
$Default = from_json($Default);
my %Default_Cpds = map { $_->{id} => 1 } @{$Default->{compounds}};

open(DB, "< ../Biochemistry/kbase_biochemistry.master.json");
my $Master="";
while(<DB>){
    $Master.=$_;
}
close(DB);
$Master = from_json($Master);

foreach my $cpd ( grep { !exists($Default_Cpds{$_->{id}}) } @{$Master->{compounds}} ){
    push(@{$Default->{compounds}},$cpd);
}

my %Test=();
foreach my $cpd (@{$Default->{compounds}}){
    print "Duplicate id: $cpd!\n" if exists($Test{$cpd->{id}});
    $Test{$cpd->{id}}=1;
}

print scalar(keys %Default_Cpds),"\t",scalar(keys %Test),"\n";

open(OUT, "> ../Biochemistry/default_biochemistry.merged_cpds.json");
print OUT to_json($Default, {pretty=>1,ascii=>1});
close(OUT);
