#!/usr/bin/perl

use warnings;
use strict;


### script to take bed12 files, and output a new bed12 file that has only the CDS for each entry


######

if (@ARGV==0) { die "\nUsage: bed12file(s).bed\n\n"; }

BEDLOOP: foreach my $bed (@ARGV) {
    if (!-e $bed) { die "\n\nterminating - cannot open bed file $bed\n\n"; }
    print "\n### working on bed file $bed\n";
    
    ## initiate output file
    my $out = $bed; $out =~ s/\.bed$//; $out .= ".justCDS.bed";
    open (OUT, "> $out");
    
    ## start getting seqs:
    open (BED, "< $bed");
    while (<BED>) {
        my $line = $_; chomp $line; my @f = split /\t/, $line;
        my $thinStart = $f[1]; my $thinEnd = $f[2];
        my $thickStart = $f[6]; my $thickEnd = $f[7];
        my $numBlocks = $f[9]; my $blockSizes = $f[10]; my $blockStarts = $f[11];
        
        ## get thickStart and thickEnd relative to the start of the whole region
        my $thickStartRelative = $thickStart - $thinStart;
        my $thickEndRelative = $thickEnd - $thinStart;
        
        ## now parse the blocks and make new blocks - just those within the CDS (thickStart-thickEnd)
        my @bStarts = split /\,/, $blockStarts; 
        my @bSizes = split /\,/, $blockSizes;
        my @newBlockStarts; my @newBlockSizes; my $newNumBlocks=0;
        for (my $i=0; $i<$numBlocks; $i++) {
            ## if we're looking for the CDS, trim to thickStart/thickEnd
            my $blockStart = $bStarts[$i];
            my $blockEnd = $bStarts[$i] + $bSizes[$i];
            my $blockSize = $bSizes[$i];
            # ignore exons entirely before thickStart
            if ($blockEnd < $thickStartRelative) { next; }
            # ignore exons entirely after thickEnd
            if ($blockStart > $thickEndRelative) { next; }
            
            # trim exons that overlap thickStart
            if ($blockStart < $thickStartRelative) {
                $blockStart = $thickStartRelative;
                $blockSize = $blockEnd - $blockStart;
            }
            # trim exons that overlap thickEnd
            if ($blockEnd > $thickEndRelative) {
                $blockEnd = $thickEndRelative;
                $blockSize = $blockEnd - $blockStart;
            }
            # block start should now be relative to thickStart rather than thinStart (as before)
            my $newBlockStart = $blockStart - $thickStartRelative;
            
            push @newBlockStarts, $newBlockStart;
            push @newBlockSizes, $blockSize;
            $newNumBlocks++;
        }
        if ($newNumBlocks == 0) {
            die "\n\nterminating - I did not find any blocks - not sure what to do. Here's the line I was processing:\n\n$line\n\n";
        }
        ## write it out
        my @newF = @f;
        # fix overall start and end
        $newF[1] = $thickStart; $newF[2] = $thickEnd;
        # fix blocks
        $newF[9] = $newNumBlocks; 
        $newF[10] = join ",", @newBlockSizes; $newF[10] .= ",";
        $newF[11] = join ",", @newBlockStarts; $newF[11] .= ",";
        my $newLine = join "\t", @newF;
        print OUT "$newLine\n";
    }
    close OUT;
    close BED;
}
