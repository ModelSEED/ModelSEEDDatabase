#!/usr/bin/env perl
use warnings;
use strict;
use JSON;
my $output;

use lib '/homes/seaver/Projects/PATRIC_Deploy/dev_container/modules/Workspace/lib/';
use lib '/homes/seaver/Projects/PATRIC_Deploy/dev_container/modules/auth/lib/';
use lib '/homes/seaver/Projects/ModelDeploy/kbapi_common/lib/';
use Bio::P3::Workspace::ScriptHelpers;
use Bio::P3::Workspace::WorkspaceClient;

my $pms_path="/chenry/public/modelsupport/biochemistry/default.biochem";
my $biochem_file="/homes/seaver/MSD_v1.0_Biochem.json";

print "Loading Biochemistry Data\n";
my $json=undef;
open(FH, "< $biochem_file");
while(<FH>){
    chomp;
    $json.=$_;
}
close(FH);
$json=from_json($json);
#print to_json($json,{pretty=>1}),"\n";

$output = Bio::P3::Workspace::ScriptHelpers::wscall("create",{ objects => [[$pms_path,"biochemistry",{},$json]], adminmode=>1, overwrite=>1 });
print "Loaded Biochemistry Data\n";
