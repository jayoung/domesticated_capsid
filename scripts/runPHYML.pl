#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

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

######## PHYML options:

my $alnType = "nt"; # nt or aa

my $subModel = "GTR";
#my $subModel = "HKY85";

## I might not want bootstraps (default for phyml is that it gives approximate Bayes branch supports). Or, perhaps I want 100 reps or more.
my $bootstrap = 0;
#my $bootstrap = "100";

my $removeBranchLengths = 1;

my $phymlOptions = "--pinv e --alpha e -f e";

my $numThreads = 1;

######## sbatch options:

### set this if on rhino and want to put out one job to each node.  Might want to switch to sbatch at some point??
my $use_sbatch = 1;

#### walltime, in days-hours (i.e. 1-0 is 1 day,   0-1 is 1 hour). default is 2 days. It might be more efficient for overall cluster usage to request something more realistic
#my $walltime = "default";
my $walltime = "0-4";

my $jobname = "phymlTree";

## get any non-default options from commandline
GetOptions("type=s" => \$alnType,
           "model=s" => \$subModel,
           "options=s" => \$phymlOptions,
           "boots=i" => \$bootstrap,
           "threads=i" => \$numThreads,
           "removeBlen=i" => \$removeBranchLengths,
           "sbatch=i" => \$use_sbatch,
           "walltime=s" => \$walltime,
           "job=s" => \$jobname
            ) or die "\n\nterminating - unknown option(s) specified on command line\n\n";


#####################
if ($use_sbatch == 1) {print "\n\nUsing sbatch to parallelize\n\n";}

if (($alnType ne "nt") & ($alnType ne "aa")) {
    die "\n\nterminating - unrecognized alignment type $alnType - should be nt or aa\n\n";
}
if (($alnType eq "aa") & ($subModel eq "GTR")) {
    die "\n\nterminating - you specified an amino acid alignment but did not change the default model from GTR, which is a nucleotide model. Try -model=JTT\n\n";
}
if (($numThreads > 1) & ($bootstrap == 0)) {
    die "\n\nterminating - you specified >1 thread but you're not bootstrapping. Not recommended\n\n";
}

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
    
    ### construct the phyml command
    if (!-e "$treeDir/$file.phyml_phyml_tree.txt") { 
        my $phymlCommand;
        my $logfile = "$file.phyml.log.txt";
        if ($numThreads == 1) {
            $phymlCommand .= "phyml ";
        } else {
            $phymlCommand .= "mpirun -n $numThreads --oversubscribe phyml-mpi ";
        }
        
        $phymlCommand .= " -i $file.phyml -d $alnType --sequential -m $subModel $phymlOptions > $logfile";
        if ($bootstrap > 0) { $phymlCommand .= " -b $bootstrap"; }
        if ($numThreads > 1) { 
            # need to double-escape if we're using sbatch:
            if ($use_sbatch == 1) { 
                $phymlCommand = "/bin/bash -c \\\"source /app/lmod/lmod/init/profile; module load OpenMPI; $phymlCommand ; module purge\\\"";
            } else {
                $phymlCommand = "/bin/bash -c \"source /app/lmod/lmod/init/profile; module load OpenMPI; $phymlCommand ; module purge\"";
            }
        }
        ## add $phymlCommand to $command
        $command .= "$phymlCommand ; ";
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

