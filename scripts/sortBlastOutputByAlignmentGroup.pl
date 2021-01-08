#!/usr/bin/perl

use warnings;
use strict;
use Cwd;

### after I did blast searches of each human Refseq, I want to sort them into groups that I will then align together.

## usage: specify the files I want to sort on the command line

#####

my $sortedDir = "sortedBlasts";
my $cwd = cwd();

##### file locations:
my $projectDir = "/fh/fast/malik_h/user/jayoung/miscMalikLab/domesticated_capsid";
my $accsTable = "getHumanPepsFromWillFiles/fullGenbankEntries/2020.23.10.humanSeqs.fullLen.Refseq.accs.txt.info.byAccession.txt";
my $groupsTable = "miscellaneous/humanGenesAlignmentGroups.txt";


###########################

if (@ARGV == 0) {
    die "\n\nTerminating - please specify some files to sort on the command line\n\n";
}

if (!-e $sortedDir) {mkdir $sortedDir;}

### read accession-to-gene table
if (!-e "$projectDir/$accsTable") {
    die "\n\nTerminating - cannot open accession-to-gene table $projectDir/$accsTable\n\n";
}
my %accs;
open (ACCS, "< $projectDir/$accsTable");
while (<ACCS>) {
    my $line = $_; chomp $line; my @f = split /\t/, $line;
    if ($f[0] eq "Accession") {next;} # header
    $accs{$f[0]} = $f[4];
}
close ACCS;

### read gene-to-alignmentGroup table
if (!-e "$projectDir/$groupsTable") {
    die "\n\nTerminating - cannot open gene-to-alignmentGroup table $projectDir/$groupsTable\n\n";
}
my %genes;
open (GENES, "< $projectDir/$groupsTable");
while (<GENES>) {
    my $line = $_; chomp $line; my @f = split /\t/, $line;
    if ($f[0] eq "Gene") {next;} # header
    $genes{$f[0]} = $f[1];
    if (!-e "$sortedDir/$f[1]") {mkdir "$sortedDir/$f[1]";}
}
close GENES;


### go through files to sort:
foreach my $file (@ARGV) {
    if (!-e $file) {
        die "\n\nTerminating - you're asking me to sort a file that does not exist: $file\n\n";
    }
    # get accession from file name
    if ($file !~ m/-/) {
        die "\n\nTerminating - script is set up to expect accession before a hyphen character but this file does not match that pattern: $file\n\n";
    }
    my $acc = (split /-/, $file)[0];
    if ($acc =~ m/\./) {  $acc = (split /\./, $acc)[0]; }
    if (!defined $accs{$acc}) {
        die "\n\nTerminating - file $file had acc $acc which is not in the accessions table\n\naccessions table: $projectDir/$accsTable\n\n";
    }
    my $gene = $accs{$acc};
    
    if (!defined $genes{$gene}) {
        die "\n\nTerminating - file $file had acc $acc gene $gene which is not in the alignmentGroup table\n\nalignmentGroup table: $projectDir/$groupsTable\n\n";
    }
    my $group = $genes{$gene};
    
    my $command = "ln -s $cwd/$file $sortedDir/$group";
    #print "Command $command\n\n";
    system($command);
}


