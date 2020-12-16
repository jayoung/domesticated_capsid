#!/usr/bin/perl

use warnings;
use strict;


### script to download chain files (for liftOver) for hg38 versus placental mammals, using the same assemblies used by the 100-way alignment


## liftOver options notes:
## these are the parameters that give the same output as the web liftOver (at least for panTro4):
# -minMatch=0.N Minimum ratio of bases that must remap. Default 0.95

## these are other parameters I toyed with, but might not use:
# -minBlocks=0.49 -fudgeThick
# -minBlocks=0.N Minimum ratio of alignment blocks or exons that must map
#              (default 1.00)
# -fudgeThick    (bed 12 or 12+ only) If thickStart/thickEnd is not mapped,
#              use the closest mapped base.  Recommended if using -minBlocks.


## I use downloaded chain files. I also tried using a remote chain file (e.g. http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToSusScr3.over.chain.gz), but it failed to find the lifted over region, despite the fact it succeeded with the same chain file stored locally. It was also slow.


### this has a list of the assemblies we expect to see
my $assemblyListFile = "/fh/fast/malik_h/grp/public_databases/UCSC/human_Dec2013/conservation_tracks/README.placentalMammals.txt";

## specify any assemblies we DON'T want to liftOver here
my %assembliesToIgnore;
$assembliesToIgnore{'hg38'}=1;  ## meaningless to chain against itself

my $chainFilesDir = "/fh/fast/malik_h/grp/public_databases/UCSC/human_Dec2013/liftOverFiles";

my $bedInputFile = "allRefSeqs.hg38.bed";

my $liftOverExecutable = "./liftOver";

my $liftOverOptions = "-minMatch=0.1";

######

## check that liftOver is in my path (need to have loaded Kent_tools module)
if (!-e $liftOverExecutable) { die "\n\nterminating - listOver executable specified in script does not exist $liftOverExecutable\n\n"; }
## old, when I was using the scicomp module:
#my $whichLiftOver = `which liftOver`;
#if ($whichLiftOver eq "") {
#    die "\n\nterminating - please run this command before running runLiftOver.pl :\n\nmodule load Kent_tools\n\n";
#}

## check input files exist
if (!-e $assemblyListFile) { die "\n\nterminating - assembly list file does not exist $assemblyListFile\n\n"; }
if (!-e $chainFilesDir) { die "\n\nterminating - chain files dir does not exist $chainFilesDir\n\n"; }
if (!-e $bedInputFile) { die "\n\nterminating - bed input file does not exist $bedInputFile\n\n"; }

open (ASSEMBLIES, "< $assemblyListFile");
my @filesGot;
while (<ASSEMBLIES>) {
    my $line = $_; chomp $line;
    print "#### assembly $line\n";
    if (defined $assembliesToIgnore{$line}) {
        print "    ignoring this one (as specified in script)\n";
        next;
    }
    my $chainFile = ucfirst($line);
    $chainFile = "hg38To".$chainFile.".over.chain.gz";
    if (!-e "$chainFilesDir/$chainFile") {
        die "\n\nterminating - cannot find chain file $chainFile\n";
    }
    
    my $outStem = $bedInputFile; 
    $outStem =~ s/\.bed$//;
    $outStem =~ s/\.hg38$//;
    
    my $outMapped = "$outStem.liftOver.$line.mapped.bed";
    my $outUnmapped = "$outStem.liftOver.$line.unmapped.bed";
    
    if (-e $outMapped) {
        print "    skipping this one, as output file exists already:$outMapped\n";
        next;
    }
    
    my $liftOverCommand = "$liftOverExecutable $liftOverOptions $bedInputFile $chainFilesDir/$chainFile $outMapped $outUnmapped";
    print "    running liftOver command:\n$liftOverCommand \n";
    system("$liftOverCommand");
}
close ASSEMBLIES;


