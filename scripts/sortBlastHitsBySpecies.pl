#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

### goes through several tab-delimited blast output files, with columns as specified in callblastjanetNRntProtquery.bioperl
# outputs a file for each species of accessions that I might get


my $maxHitLength = 10000; ## to eliminate genomic seqs
my $minQueryLengthCoverage = 0.5; ## e.g.must match >50% of length of query
my $minIdent = 0.4; ## 40% amino acid identity
my $outDir = "hits_each_species";

## get any non-default options from commandline
GetOptions("maxHlen=i" => \$maxHitLength,
           "minQcov=i" => \$minQueryLengthCoverage,
           "minIdent=i" => \$minIdent,
           "outDir=s" => \$outDir
            ) or die "\n\nterminating - unknown option(s) specified on command line\n\n";

#####

### collect all the hits, by species
my %hitBySpecies;
foreach my $file (@ARGV){
    if (!-e $file) { die "\n\nterminating - cannot open file $file\n\n"; }
    print "\n####### working on file $file\n";
    open (IN, "< $file");
    my $firstLine = "yes";
    while(<IN>) {
        my $line = $_; chomp $line;
        ## check that header looks correct
        if ($firstLine eq "yes") {
            $line =~ s/\t/ /g;
            ## this comes from callblastjanetNRntProtquery.bioperl
            my $shouldBe = "qaccver qlen saccver slen evalue score bitscore pident ppos length gapopen qstart qend sstart send sframe staxid ssciname scomname stitle";
            if ($line ne $shouldBe) {
                die "\n\nterminating - header does not look like I expect it to - was this generated using the callblastjanetNRntProtquery.bioperl script?\n\n";
            }
            $firstLine = "no";
            next;
        }
        ### other lines
        my @f = split /\t/, $line;
        ## filtering
        my $hitLen = $f[3]; if ($hitLen > $maxHitLength) {next;}
        my $ident = $f[7]; if ($ident < $minIdent) {next;}
        my $queryCov = $f[9] / $f[1]; if ($queryCov < $minQueryLengthCoverage) {next;}
        ## keep the rest
        my $species = $f[17]; $species =~ s/\s/_/g;
        my $acc = $f[2];
        $hitBySpecies{$species}{$acc}=1;
    }
    close IN;
}

### output hits for each species
if (!-e $outDir) {mkdir $outDir;}
foreach my $species (sort keys %hitBySpecies) {
    my $out = "$outDir/hits_$species.txt";
    my @accs = sort keys %{$hitBySpecies{$species}};
    open (OUT, "> $out");
    foreach my $acc (@accs) { print OUT "$acc\n"; }
    close OUT;
}
