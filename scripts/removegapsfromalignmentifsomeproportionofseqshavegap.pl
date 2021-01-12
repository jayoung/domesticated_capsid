#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;

###removes columns from alignment if all sequences have a gap at that position


## do we treat Ns as gaps?
my $Nisgap = 0;

### remove position if this proportion of sequences have a gap
my $proportion = 0.5;

## get any non-default options from commandline
GetOptions("threshold=s" => \$proportion,
           "Ngap=i" => \$Nisgap
            ) or die "\n\nterminating - unknown option(s) specified on command line\n\n";

########################################

if ($Nisgap == 1) {
   print "\n\ncounting Ns as gaps - if this is a protein alignment that's wrong\n\n";
} else {print "Ns are not counted as gaps\n";}

print "removing columns from the alignment if more than $proportion of them have gaps\n";

foreach my $infile (@ARGV) {
    if (!-e $infile) { die "\n\nterminating - file $infile not found\n\n";}
    #my $outfile = "$infile" . ".severelydegapped_" . $proportion;
    my $outfile = $infile; $outfile =~ s/\.fasta$//; $outfile =~ s/\.fa$//;
    $outfile .= ".degap_" . $proportion . ".fa";
    print "processing file $infile\n";
    my $seqIN = Bio::SeqIO->new(-file=>"< $infile", -format=>"fasta");
    my @seqs;
    while (my $seq = $seqIN->next_seq) { push @seqs, $seq; }
    my $numseqsneedgaptoremove = @seqs;
    $numseqsneedgaptoremove = $numseqsneedgaptoremove * $proportion;
    
    my @seqorder;
    #go through sequences (the hash sequences will be a hash of arrays, one letter per position)
    my $seqOUT = Bio::SeqIO->new(-file=>">$outfile",-format=>'fasta');
    my $length;
    my %sequences;
    foreach my $seq (@seqs) {
        my $seqname = $seq->display_id;
        my $sequence = $seq->seq;
        push @seqorder, $seqname;
        $length = length ($sequence);
        for (my $i=0;$i<=$length;$i++) {
            push @{$sequences{$seqname}}, substr($sequence,$i,1);
        }
    }
    
    #now go through each position in alignment and if at least one non-gap character, add to string in the new hash, newseqs
    my %newseqs;
    foreach my $key1 (keys %sequences) { $newseqs{$key1} = ""; }
    for (my $j=0;$j<=$length;$j++) {
        my $countgap = 0;
        foreach my $key2 (keys %sequences) {
            if (@{$sequences{$key2}}[$j] eq "-") {$countgap++;}
            if ((@{$sequences{$key2}}[$j] eq "N") & ($Nisgap == 1)) {$countgap++;}
        }
        if ($countgap <= $numseqsneedgaptoremove ) {
            foreach my $key3 (keys %sequences) {
                $newseqs{$key3} .= @{$sequences{$key3}}[$j];
            } 
        }
    }
    #### write output
    foreach my $key4 (@seqorder) {
        my $seqObj = Bio::Seq->new(-seq=>"$newseqs{$key4}", -display_id=>"$key4");
        $seqOUT -> write_seq($seqObj);
    }
}
