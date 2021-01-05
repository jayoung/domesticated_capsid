#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

## xx to do:
# capture slurm output in a log file
# incorporate phyml-mpi option

####### runPHYMLonNucleotideAlignmentAndRemoveBranchLengths.pl alignmentFiles.fa

#### e.g. for the file H2Aaln1.fa.trim49_294
# mkdir H2Aaln1.fa.trim49_294_PHYMLtree
# cp H2Aaln1.fa.trim49_294 H2Aaln1.fa.trim49_294_PHYMLtree
# cd H2Aaln1.fa.trim49_294_PHYMLtree
# fasta2phyml.pl H2Aaln1.fa.trim49_294 
# phyml -i H2Aaln1.fa.trim49_294.phyml -d nt --sequential -m GTR --pinv e --alpha e -f e
# changenamesinphyliptreefileguessaliasfilename.pl   H2Aaln1.fa.trim49_294.phyml_phyml_tree
# removebranchlengthsAndBootstrapsFromphyliptree.pl H2Aaln1.fa.trim49_294.phyml_phyml_tree.names
### final outfile will be called H2Aaln1.fa.trim49_294.phyml_phyml_tree.names.nolen

##### I might not want bootstraps (default for phyml is that it gives approximate Bayes branch supports). Or, perhaps I want 100 reps or more.
my $bootstrap = 0;
#my $bootstrap = "100";

my $subModel = "GTR";
#my $subModel = "HKY85";

my $removeBranchLengths = 1;


### set this if on rhino and want to put out one job to each node.  Might want to switch to sbatch at some point??
my $use_sbatch = 1;

#### walltime, in days-hours (i.e. 1-0 is 1 day,   0-1 is 1 hour). default is 2 days. It might be more efficient for overall cluster usage to request something more realistic
#my $walltime = "default";
my $walltime = "0-4";

my $jobname = "phymlTree";

## get any non-default options from commandline
GetOptions("model=s" => \$subModel,
           "boots=i" => \$bootstrap,
           "removeBlen=i" => \$removeBranchLengths,
           "sbatch=i" => \$use_sbatch,
           "walltime=s" => \$walltime,
           "job=s" => \$jobname
            ) or die "\n\nterminating - unknown option(s) specified on command line\n\n";


#####################
if ($use_sbatch == 1) {print "\n\nUsing sbatch to parallelize\n\n";}

foreach my $file (@ARGV) {
    if (!-e $file) {
        die "\n\nterminating - file $file does not exist\n\n";
    }
    my $treeDir = "$file"."_PHYMLtree";
    my $outfile = "$treeDir/$file.phyml_phyml_tree.names.nolen";
    if (-e $outfile) {
        print "\nskipping $file - outfile $outfile exists already\n\n";
        next;
    }
    my $command = "";
    # mkdir H2Aaln1.fa.trim49_294_PHYMLtree
    if (!-e $treeDir) { $command .= "mkdir $treeDir ; ";}
    # cp H2Aaln1.fa.trim49_294 H2Aaln1.fa.trim49_294_PHYMLtree
    if (!-e "$treeDir/$file") { $command .= "cp $file $treeDir ; ";}
    # cd H2Aaln1.fa.trim49_294_PHYMLtree
    $command .= "cd $treeDir ; ";
    # fasta2phyml.pl H2Aaln1.fa.trim49_294 
    if (!-e "$treeDir/$file.phyml") { $command .= "fasta2phyml.pl $file ; "; }
    # phyml -i H2Aaln1.fa.trim49_294.phyml -d nt --sequential -m $subModel --pinv e --alpha e -f e
    if (!-e "$treeDir/$file.phyml_phyml_tree.txt") { 
        my $logfile = "$file.phyml.log.txt";
        $command .= "phyml -i $file.phyml -d nt --sequential -m $subModel --pinv e --alpha e -f e > $logfile";
        if ($bootstrap > 0) { $command .= " -b $bootstrap"; }
        $command .= " ; "; 
    }
    $command .= "changenamesinphyliptreefileguessaliasfilename.pl $file.phyml_phyml_tree.txt ; ";
    if ($removeBranchLengths == 1) {
        $command .= "removebranchlengthsAndBootstrapsFromphyliptree.pl $file.phyml_phyml_tree.txt.names";
    }
    
    if ($use_sbatch == 1) {
        my $time = "";
        if ($walltime ne "default") { $time = "-t $walltime"; }
        $command = "sbatch --job-name=$jobname $time --wrap=\"$command\"";
    }
    print "\n\ncommand $command\n\n";
    system($command);
}

if ($use_sbatch == 1) {
    print "\n\nSet all jobs going - use sq command to monitor whether there are still any phymlTree commands running\n\n";
}

