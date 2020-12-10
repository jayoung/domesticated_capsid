#!/usr/bin/perl

use warnings;
use strict;

#### this works on any set of seqs. Uses default scoring system (frameshifts and stops quite penalized)

### set this if on rhino and want to put out one job to each node
my $use_sbatch = "yes";
#my $use_sbatch = "no";

#### days-hours (i.e. 1-0 is 1 day,   0-1 is 1 hour). default is 2 days. It might be more efficient for overall cluster usage to request something more realistic
#my $walltime = "default";
my $walltime = "0-12";

#### if set to yes, it'll just print the command up on the screen - won't actually do it.
my $debug = "no";
#my $debug = "yes";


##################

if ($use_sbatch eq "yes") {print "\n\nUsing sbatch to parallelize\n\n";}

foreach my $file (@ARGV) {
    if (!-e $file) { die "\n\nterminating - file $file does not exist\n\n";}
    print "\n#### working on file $file\n";
    my $hits = $file;  
    my $out = $file;
    $out =~ s/\.fa$//;
    $out .= "_aln1";
    my $outNT = "$out"."_NT.fa";
    my $outAA = "$out"."_AA.fa";
    my $outLog = "$out"."_MACSE.log.txt";
    if (-e $outNT) {
        print "    Skipping file $file - outfile $outNT exists already\n\n";
        next;
    }
    
    my $command = "java -jar -Xmx600m /fh/fast/malik_h/grp/malik_lab_shared/bin/macse_v1.2.jar -prog alignSequences -seq $file -out_NT $outNT -out_AA $outAA > $outLog";
    
    if ($use_sbatch eq "yes") {
        my $time = "";
        if ($walltime ne "default") { $time = "-t $walltime"; }
        $command = "/bin/bash -c \\\"source /app/lmod/lmod/init/profile; module load Java/1.8.0_181 ; $command\\\"";
        $command = "sbatch $time --job-name=MACSE --wrap=\"$command\"";
    } 
    print "command:\n    $command\n\n";
    if ($debug eq "no") { system($command); }
}

if ($use_sbatch eq "yes") {
    print "\n\nSet all jobs going - use sq command to monitor whether there are still any MACSE commands running\n\n";
}