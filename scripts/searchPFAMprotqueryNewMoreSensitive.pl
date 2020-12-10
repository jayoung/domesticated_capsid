#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

### runs pfamscan_lwp.pl 

#### still some work to do to make this a nice tidy script:

### does not do any checking for whether we've done this search already
### does not try to combine output files

my $email = "jayoung\@fhcrc.org";

foreach my $file (@ARGV) {
    if (!-e $file) {
        die "\n\nterminating - cannot open file $file\n\n";
    }
    
    my $outDir = $file . "_searchPFAM_moreSensitive";
    
    ## go through the input file to see what query seqs were present - check there are no duplicate seqnames
    my $seqIN = Bio::SeqIO->new(-file=>"< $file", -format=>"fasta");
    my %seqsSeen;
    my $numSeqsWithoutOutputAlready = 0;
    while(my $seq=$seqIN->next_seq) {
        #print "seq $seq\n"; die;
        my $seqname = $seq->display_id();
        if (defined $seqsSeen{$seqname}) {
            die "\n\nterminating  - problem: the seq file contains two seqs with the same name. pfamscan_lwp.pl will overwrite results if you run it with these queries\n\n";
        }
        $seqsSeen{$seqname}=1;
        $numSeqsWithoutOutputAlready++;
        my $tempoutfile = $seqname;
        $tempoutfile =~ s/\./_/g;
        $tempoutfile =~ s/-/_/g;
        $tempoutfile =~ s/[\(\)]/_/g;
        $tempoutfile = "$outDir/$tempoutfile.out.txt"; 
        #print "checking for out file $tempoutfile blah\n";
        if (-e $tempoutfile) {
            $numSeqsWithoutOutputAlready--;
        } else {
            print "    did not find outfile $tempoutfile\n";
        }
    }
    #print "numSeqsWithoutOutputAlready $numSeqsWithoutOutputAlready\n"; die;
    ## make a dir wheere the results will go and cd into it
    if (!-e $outDir) {mkdir $outDir;}
    chdir $outDir;
    ## make a link within that dir to the query file:
    if (!-e $file) { system("ln -s ../$file ."); }
    
    ## run the search
    ## the result will be two files per query seq that was in the query file:  
    # (a) a .out.txt - I'll capture this output in a combined output file
    # (b) a .sequence.txt file (simply the individual query seq) - I'll delete this
    if ($numSeqsWithoutOutputAlready > 0) {
        my $command = "pfamscan_lwp.pl --evalue 10 --multifasta --useSeqId --email $email --format txt $file";
        print "\n\n##### input file $file - running command $command\n";
        system($command);
    } else {
        print "\n\n##### input file $file - output files already exist for all queries, concatenating existing results\n\n";
    }
    ## capture the output in a single combined file
    my $outFile = $file . ".PFAMscan.txt";
    open (OUT, "> $outFile");
    my $printedHeaderYet = "no";
    foreach my $seq (sort keys %seqsSeen) {
        my $seqOutFile = $seq; 
        $seqOutFile =~ s/\./_/g;
        $seqOutFile =~ s/-/_/g;
        $seqOutFile =~ s/[\(\)]/_/g;
        $seqOutFile .= ".sequence.txt";
        if (-e $seqOutFile) { unlink $seqOutFile; }
        my $PFAMoutFile = $seq; 
        $PFAMoutFile =~ s/\./_/g;
        $PFAMoutFile =~ s/-/_/g;
        $PFAMoutFile =~ s/[\(\)]/_/g;
        $PFAMoutFile .= ".out.txt";
        if (!-e $PFAMoutFile) {
            die "\n\nterminating - expected an output file called $PFAMoutFile - something is wrong\n\n";
        }
        open (IN, "< $PFAMoutFile");
        while (<IN>) {
            my $line = $_; chomp $line;
            if (($line =~ m/^#/) & ($line !~ m/^#\s<seq\sid>/)) { next; }
            if ($line eq "") {next;}
            if (($line =~ m/^#\s<seq\sid>/) & ($printedHeaderYet eq "yes")) { next; }
            if (($line =~ m/^#\s<seq\sid>/) & ($printedHeaderYet eq "no")) {
                $line =~ s/^#\s//;
                $line =~ s/>\s</\t/g;
                $line =~ s/>$//;
                $line =~ s/^<//;
                print OUT "$line\n";
                $printedHeaderYet = "yes";
                next;
            }
            $line =~ s/\s+/\t/g;
            $line =~ s/\t$//;
            print OUT "$line\n";
        }
        close IN;
    }
    close OUT;
    chdir "..";
}
