#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header=1;

#Start with list of formulas with differing protons
open(FH, "< ../Comparisons/Formula_Differences.txt");
$header=1;
my %Cpds=();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next unless $temp[0] eq "NUMBERS";
    next unless $temp[2] =~ /^H:\d+$/;
    $Cpds{$temp[1]}=$_;
}
close(FH);

#Add aliases for KEGG and MetaCyc to match InChI structures
open(FH, "< ../Aliases/KEGG.aliases");
my %Aliases = ();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next unless exists($Cpds{$temp[2]});

    $Aliases{$temp[0]} = (split(/\|/,$temp[2]))[0];
}
close(FH);

open(FH, "< ../Aliases/MetaCyc.aliases");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next unless exists($Cpds{$temp[2]});

    $Aliases{$temp[0]} = (split(/\|/,$temp[2]))[0];
}
close(FH);

use InChIs;

#Collect structures
open(FH, "< ../Structures/KEGG_Charged_InChI.txt");
my %InChIs = ();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next unless exists($Aliases{$temp[0]});

    my ($formula,$layers)=InChIs::parse($temp[1]);
    $InChIs{$Aliases{$temp[0]}}{OLD}=$formula;
    $formula = InChIs::adjust_protons($formula,$layers->{p});
    $InChIs{$Aliases{$temp[0]}}{NEW}=$formula;

}
close(FH);

open(FH, "< ../Structures/MetaCyc_Charged_InChI.txt");
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    next unless exists($Aliases{$temp[0]}) && !exists($InChIs{$Aliases{$temp[0]}});

    my ($formula,$layers)=InChIs::parse($temp[1]);
    $InChIs{$Aliases{$temp[0]}}{OLD}=$formula;
    $formula = InChIs::adjust_protons($formula,$layers->{p});
    $InChIs{$Aliases{$temp[0]}}{NEW}=$formula;
}
close(FH);
    
use Bio::KBase::fbaModelServices::ScriptHelpers qw( getToken );
my $AToken = getToken();

use Bio::KBase::fbaModelServices::Impl;
my $FBAImpl = Bio::KBase::fbaModelServices::Impl->new({'fbajobcache' => "/homes/seaver/Projects/KBase_Scripts/FBA_Scripts/JobsCache",
						       'jobserver-url' => "http://kbase.us/services/workspace",
						       'fbajobdir' => "/tmp/fbajobs",
						       'mfatoolkitbin' => "/vol/model-prod/kbase/MFAToolkit/bin/mfatoolkit",
#						       'mfatoolkitbin' => "/homes/seaver/Software/MFAToolkit/bin/mfatoolkit",
						       'probanno-url' => "http://140.221.85.86:7073/",
						       'mssserver-url' => "http://bio-data-1.mcs.anl.gov/services/ms_fba",
						       'accounttype' => "kbase",
						       'workspace-url' => "http://kbase.us/services/ws",
						       'defaultJobState' => "queued",
						       'gaserver-url' => "http://kbase.us/services/genome_annotation",
						       'idserver-url' => "http://kbase.us/services/idserver"});
$FBAImpl->_setContext(undef,{auth=>$AToken});

my $bioObj = $FBAImpl->_get_msobject("Biochemistry","kbase","plantdefault");
foreach my $cpd (sort keys %InChIs){
    my $obj = $bioObj->getObject("compounds",$cpd);
    print $cpd,"\t",$obj->formula(),"\t",$InChIs{$cpd}{OLD},"\t",$InChIs{$cpd}{NEW},"\n" if $obj->formula() eq $InChIs{$cpd}{NEW};
}
	
