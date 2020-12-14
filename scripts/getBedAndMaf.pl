#!/usr/bin/perl

use warnings;
use strict;


## gets refseq accessions from a table, gets a single line bed12 file giving coordinates of that gene in hg38, and runs mafFrags to get the maf for just placental mammals for each gene.  Run without arguments - all input files are specified below:

# if there's >1 entry for this RefSeq ID in the whole-genome bed file, I try to get rid of alt chromosomes etc

my $baseDir = "/fh/fast/malik_h/user/jayoung/miscMalikLab/domesticated_capsid";

my $accsFile = "getHumanPepsFromWillFiles/fullGenbankEntries/2020.23.10.humanSeqs.fullLen.Refseq.accs.txt.info.byAccession.txt";

my $wholeGenomeBedFile = "/home/jayoung/public_databases/UCSC/human_Dec2013/misc_tracks/hg38.ncbiRefSeq.2020dec3.bed";

my $speciesListFile = "/fh/fast/malik_h/grp/public_databases/UCSC/human_Dec2013/conservation_tracks/README.placentalMammals.txt";

my $outFileSuffix = ".hg38.100way_placMamm.maf";

###############

### check all input files exist
if (!-e "$baseDir") {
    die "\n\nterminating - baseDir does not exist $baseDir\n\n";
}
if (!-e "$baseDir/$accsFile") {
    die "\n\nterminating - accsFile does not exist $baseDir/$accsFile\n\n";
}
if (!-e $wholeGenomeBedFile) {
    die "\n\nterminating - wholeGenomeBedFile does not exist $wholeGenomeBedFile\n\n";
}
if (!-e $speciesListFile) {
    die "\n\nterminating - speciesListFile does not exist $speciesListFile\n\n";
}

# check that mafFrags is in my path (need to have loaded Kent_tools module)
my $whichMafFrags = `which mafFrags`;
if ($whichMafFrags eq "") {
    die "\n\nterminating - please run this command before running getBedAndMaf.pl :\n\nmodule load Kent_tools\n\n";
}

# work on each accession
open (ACCS, "< $baseDir/$accsFile");
while (<ACCS>) {
    my $line = $_; 
    if ($line =~ m/^Accession/) { next; }
    my @f = split /\t/, $line;
    my $acc = $f[3]; my $gene = $f[4]; 
    print "\n## gene $gene acc $acc\n";
    my $bedFile = $gene . "_" . $acc . ".bed";
    my $mafFile = $gene . "_" . $acc . $outFileSuffix;
    if (-e $mafFile) {
        print "    skipping this one - maf file exists already\n"; 
        next;
    }
    
    # grep to get a bed file for a single transcript
    my $grepCommand = "grep \'$acc\' $wholeGenomeBedFile > $bedFile";
    `$grepCommand`;
    # test for empty output of grep and report
    if (!-e $bedFile) { die "\n\nterminating - grep command failed\n\n"; }
    if (-z $bedFile) { die "\n\nterminating - grep command gave empty output\n\n"; }
    
    ### check whether there was >1 bed entry, if so, do something about it
    open (BED, "< $bedFile");
    my @lines = <BED>;
    close BED;
    my $numLines = @lines;
    if ($numLines > 1) {
        print "    WARNING - found >1 bed file entry, taking only one of them\n";
        
        ## choose one entry, hopefully one that's not on an alt/random/Un chromosome
        my @goodLines;
        my $lineToKeep;
        foreach my $line (@lines) {
            my $chr = (split /\t/, $line)[0];
            if ($chr =~ m/_alt|_random|chrUn/) {next;}
            push @goodLines, $line;
        }
        my $numGoodLines = @goodLines;
        if ($numGoodLines == 0) {
            print "        WARNING 2a - found no entries on a good chromosome\n";
            $lineToKeep = $lines[0];
        }
        if ($numGoodLines > 1) {
            print "        WARNING 2b - found >1 entry on a good chromosome\n";
            $lineToKeep = $goodLines[0];
        }
        if ($numGoodLines == 1) {
            print "        successfully filtered to just the entry on a good chromosome\n";
            $lineToKeep = $goodLines[0]; 
        }
        open (BED, "> $bedFile");
        print BED $lineToKeep;
        close BED;
    }
    
    # now run mafFrags
    my $mafFragsCommand = "mafFrags -orgs=$speciesListFile -bed12 -thickOnly hg38 multiz100way $bedFile $mafFile";
    print "running mafFrags command:\n\n$mafFragsCommand\n";
    `$mafFragsCommand`;
    # test for empty output of mafFrags and report
    if (!-e $mafFile) { die "\n\nterminating - mafFrags command failed\n\n"; }
    if (-z $mafFile) { die "\n\nterminating - mafFrags command gave empty output\n\n"; }
}
close ACCS;

