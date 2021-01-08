#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;
use Cwd;

###### groups sequences and aligns them

# input = ungapped seqs that might be orthologs to what we're studying
# blast each on back to a database of human refseq mRNAs
# record the best blast hit
# uses a lookup table to decide which alignment group it goes in
# make an output file for each alignment group 
# do the alignment


##### file locations so I can get gene names and groups:
my $projectDir = "/fh/fast/malik_h/user/jayoung/miscMalikLab/domesticated_capsid";
my $groupsTable = "miscellaneous/humanGenesAlignmentSubgroups.txt";
my $callMACSEscript = "scripts/runMACSEsimpleInFrameAlignment.pl";


###########################

my $refSeqDB = "/fh/fast/malik_h/grp/public_databases/UCSC/human_Dec2013/misc_tracks/refMrna.2020dec23.geneNames.fa";

my $baseDir = "/fh/fast/malik_h/user/jayoung/miscMalikLab/domesticated_capsid";

my $numThreads = 4;

##################

### read gene-to-alignmentGroup table
if (!-e "$projectDir/$groupsTable") {
    die "\n\nTerminating - cannot open gene-to-alignmentGroup table $projectDir/$groupsTable\n\n";
}
my %geneGroups;
open (GENES, "< $projectDir/$groupsTable");
while (<GENES>) {
    my $line = $_; chomp $line; my @f = split /\t/, $line;
    if ($f[0] eq "Gene") {next;} # header
    $geneGroups{$f[0]} = $f[1];
}
close GENES;


### get a list of seqs in $refSeqDB
if (!-e $refSeqDB) {
    die "\n\nterminating - refSeqDB file $refSeqDB does not exist\n\n";
}
print "\nGetting names of seqs in database\n";
my $db = Bio::SeqIO->new(-file=>"< $refSeqDB", -format=>"fasta");
my %dbSeqs;
while (my $db_seq = $db->next_seq()) {
    #my $id = $db_seq->display_id();
    my $id = $db_seq->description();
    $dbSeqs{$id} = 1;
}
my $num = keys %dbSeqs;
print "    read $num database seqs\n";

my @warnings; 

foreach my $file (@ARGV) {
    if (!-e $file) {
        die "\n\nterminating - file $file does not exist\n\n";
    }
    print "####### working on file $file\n";
    ## make a subdir to work in, cd into it
    my $dir = $file . "_blastBackCheck";
    if (!-e $dir) {mkdir $dir;}
    chdir $dir;
    if (!-e $file) { system("cp ../$file ."); }
        
    ## get ungapped seqs (not really necessary - my input is probably ungapped, but I'll keep it here as I'm too lazy, and maybe I'll use gapped input at some point)
    my %seqs;  ## hash - key = seqname, value= bioperl seq object
    my %gappedSeqs;  ## hash - key = seqname, value= bioperl seq object
    my $seqIN = Bio::SeqIO->new(-file=>"< $file", -format=>"fasta");
    while (my $seq = $seqIN->next_seq()) {
        my $id = $seq->display_id();
        ## remove gaps
        my $letters = $seq->seq();
        my $gappedLetters = $letters;  $gappedLetters =~ s/!/-/g;
        ## make a new object containing the gapped seq and store it.  A plain old copy did not work (both objects got updated). I need to do a "deep copy" instead - I think I could make it work with code like this:
        # use Storable qw(dclone);
        # $r2 = dclone($r1); 
        # https://docstore.mik.ua/orelly/perl4/cook/ch11_13.htm
        my $desc = $seq->description();
        #print "seq $id desc $desc\n";
        my $gappedSeq = Bio::Seq->new(-display_id=>$id, -seq=>$gappedLetters, -description=>$desc);
        $gappedSeqs{$id} = $gappedSeq;
        
        $letters =~ s/-//g; $letters =~ s/!//g;
        $seq->seq($letters);
        ## store seq
        $seqs{$id} = $seq;
    }
    
    my %seqsByGroup;
    my $reportFile = "$file.blastCheck.txt";
    my $reportFile2 = "$file.blastCheckSummary.txt";
    open (OUT, "> $reportFile");
    print OUT "Seq\tGroup\tBestBlastHitGeneName\tBestBlastHitAcc\n";
    
    ## go through seqs, blast to $refSeqDB
    my $blastDir = "blastBacks";
    if (!-e $blastDir) { mkdir $blastDir; }
    my $matchCount = 0; my $noMatchCount = 0;
    foreach my $querySeq (sort keys %seqs) {
        ## write seq file
        my $queryFile = "$blastDir/$querySeq.fa";
        my $seqOUT = Bio::SeqIO->new(-file=>"> $queryFile", -format=>"fasta");
        $seqOUT->write_seq($seqs{$querySeq});
        ## do blast
        my $blastOutFile = "$queryFile.blastnRefDB";
        my $blastCommand = "blastn -task blastn -db $refSeqDB -query $queryFile -out $blastOutFile -num_alignments 50 -num_descriptions 50 -evalue 0.00001 -num_threads $numThreads";
        if (!-e $blastOutFile) {
            print "    Running blast on $querySeq\n";
            #print "        $blastCommand\n";
            system($blastCommand);
        }
        ## parse blast
        my $parsedBlast = "$blastOutFile.procnew.txt";
        if (!-e $parsedBlast) {
            my $parseCommand = "$baseDir/scripts/blastparsenew.bioperl $blastOutFile";
            #print "    Parsing blast\n";
            system($parseCommand);
        }
        
        ## look at parsed blast. our human seq should be top hit or have equal score to top hit. output will be a table with one row per seq in the alignment, columns are seqname, then yes or no about blast back (and if no, what other seq does it blast back to) 

        my @equallyBestHits;
        open (IN, "< $parsedBlast");
        my $matchedGeneName;
        my $matchedGeneAcc;
        my $geneGroup = "unknown";
        while (<IN>) {
            my $line = $_; chomp $line;
            if ($line =~ m/^QUERY/) {next;}
            if ($line =~ m/^-----/) {next;}
            my @f = split /\t/, $line;
            $matchedGeneAcc = $f[2];
            if ($f[21] ne "") {
                $matchedGeneName = $f[21];
            } else {
                $matchedGeneName = $matchedGeneAcc;
            }
            if (defined $geneGroups{$matchedGeneName}) {
                $geneGroup = $geneGroups{$matchedGeneName}; 
            }
            last; ## I only use the first hit
        }
        close IN;
        push @{$seqsByGroup{$geneGroup}}, $seqs{$querySeq};
        print OUT "$querySeq\t$geneGroup\t$matchedGeneName\t$matchedGeneAcc\n"; 
    } # end of going through query seqs
    close OUT;
    
    open (OUT2, "> $reportFile2");
    my $numSeqs = keys %seqs;
    print "\n### done. There were $numSeqs seqs in the file. Of those:\n";
    print OUT2 "\n### done. There were $numSeqs seqs in the file. Of those:\n";
    
    ### go through each group, make output file
    print "\nGroup\tNumSeqs\n";
    print OUT2 "Group\tNumSeqs\n";
    foreach my $group (sort keys %seqsByGroup) {
        my @seqs = @{$seqsByGroup{$group}};
        my $num = @seqs;
        #print "$group\t$num\n";
        print OUT2 "$group\t$num\n";
        my $groupSeqsFile = $file;
        $groupSeqsFile =~ s/\.fasta$//; $groupSeqsFile =~ s/\.fa$//;
        $groupSeqsFile .= ".group_$group.fa";
        my $groupSeqsOUT = Bio::SeqIO->new(-file=>"> $groupSeqsFile", -format=>"fasta");
        foreach my $seq (@seqs) { $groupSeqsOUT->write_seq($seq);}
        
        # submit MACSE job 
        if ($group ne "unknown") {
            my $command = "$projectDir/$callMACSEscript $groupSeqsFile";
            system($command);
        }
    }
    print "\n\n";
    system("cat $reportFile2");
    print "\n\n";
    ## cd back to the starting dir
    chdir "..";
    close OUT2;
}


