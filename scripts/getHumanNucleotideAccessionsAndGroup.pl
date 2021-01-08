#!/usr/bin/perl

use warnings;
use strict;
use Cwd;

### to start off my alignments, I want the human nucleotide refseq accessions that correspond to the protein accessions that Will got, and I want to group them according to which alignment I'll put them in


##### file locations:
my $projectDir = "/fh/fast/malik_h/user/jayoung/miscMalikLab/domesticated_capsid";
my $accsTable = "getHumanPepsFromWillFiles/fullGenbankEntries/2020.23.10.humanSeqs.fullLen.Refseq.accs.txt.info.byAccession.txt";
my $groupsTable = "miscellaneous/humanGenesAlignmentGroups.txt";


###########################

my $cwd = cwd();
if ($cwd ne "$projectDir/alignments") {
    die "\n\nTerminating - please run this from the $projectDir/alignments dir. You are in $cwd\n\n";
}



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
}
close GENES;

### read accession-to-gene table, now looking at nucleotide accs
if (!-e "$projectDir/$accsTable") {
    die "\n\nTerminating - cannot open accession-to-gene table $projectDir/$accsTable\n\n";
}
my %accsByGroup;
open (ACCS, "< $projectDir/$accsTable");
while (<ACCS>) {
    my $line = $_; chomp $line; my @f = split /\t/, $line;
    if ($f[0] eq "Accession") {next;} # header
    my $nucAcc = $f[3];
    my $gene = $f[4];
    if (!defined $genes{$gene}) {
        die "\n\nTerminating - gene $gene is not in a group\n\n";
    }
    my $group = $genes{$gene};
    push @{$accsByGroup{$group}}, $nucAcc;
}
close ACCS;

### output accessions to get for each group
foreach my $group (sort keys %accsByGroup) {
    if ($group eq "none") {next;}
    print "Working on group $group\n";
    if (!-e $group) {mkdir $group;}
    my $accsFile = "$group/$group" . "_accs.txt";
    open (OUT, "> $accsFile");
    foreach my $acc (@{$accsByGroup{$group}}) { print OUT "$acc\n"; }
    close OUT;
}

