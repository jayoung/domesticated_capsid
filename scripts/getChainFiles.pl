#!/usr/bin/perl

use warnings;
use strict;


### script to download chain files (for liftOver) for hg38 versus placental mammals, using the same assemblies used by the 100-way alignment

### this has a list of the assemblies we expect to see
my $assemblyListFile = "/fh/fast/malik_h/grp/public_databases/UCSC/human_Dec2013/conservation_tracks/README.placentalMammals.txt";

## specify any assemblies we DON'T want chain files for here
my %assembliesToIgnore;
$assembliesToIgnore{'hg38'}=1;  ## meaningless to chain against itself

my $destinationDir = "/fh/fast/malik_h/grp/public_databases/UCSC/human_Dec2013/liftOverFiles";



######

if (!-e $assemblyListFile) { die "\n\nterminating - assembly list file does not exist $assemblyListFile\n\n"; }
if (!-e $destinationDir) { die "\n\nterminating - destination dir does not exist $destinationDir\n\n"; }

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
    if (-e "$destinationDir/$chainFile") {
        print "    already got chain file for this one\n";
        next;
    }
    my $ftpAddress = "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/$chainFile";
    my $downloadCommand = "wget --directory-prefix=$destinationDir \'$ftpAddress\'";
    print "    running command\n$downloadCommand\n\n";
    system("$downloadCommand");
    push @filesGot, $chainFile;
}
close ASSEMBLIES;

## now check all files succeeded
print "\n\n######## Done. Checking files\n\n";
my $numProblems = 0;
foreach my $file (@filesGot) {
    if (!-e "$destinationDir/$file") { print "    FAILED! $file\n"; $numProblems++; }
    if (-z "$destinationDir/$file") { print "    EMPTY! $file\n"; $numProblems++; }
}
if ($numProblems==0) {
    print "\n    Got all chain files with no problems\n\n";
}
