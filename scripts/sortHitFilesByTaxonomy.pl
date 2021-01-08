#!/usr/bin/perl

use warnings;
use strict;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;

## takes files with names like hits_Sus_scrofa.txt, and a grouping file, and puts files into separate directories

my $groupingsFile = "/fh/fast/malik_h/user/jayoung/miscMalikLab/domesticated_capsid/miscellaneous/taxonomicGroupsDesired.txt";

########################### script starts here - do not edit below this line

if (@ARGV == 0) {
    die "\n\nterminating - please specify hits*txt files\n\n";
}
if (!-e $groupingsFile) {
    die "\n\nterminating - cannot find the groupings file $groupingsFile\n\n";
}


## open database connection using local flatfile.
# Directory specifies where the indexes will go, or where they will be found. Takes 2 mins to run if indices do not exist, but then is very fast after that
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile',
      -nodesfile => '/fh/fast/malik_h/grp/public_databases/NCBI/taxonomy/nodes.dmp',
      -namesfile => '/fh/fast/malik_h/grp/public_databases/NCBI/taxonomy/names.dmp',
      -directory => '/fh/fast/malik_h/grp/public_databases/NCBI/taxonomy');
my $tree_functions = Bio::Tree::Tree->new();


print "Reading groupings file\n";
open (GROUPS, "< $groupingsFile");
my %groupings;
my $anyWarnings = 0;
while (<GROUPS>) {
    my $line = $_; chomp $line; my @f = split /\t/, $line;
    if ($f[0] eq "Group my name") {next;} # header
    #print "#### line $line\n";
    ## check include exists in taxonomy database
    my $include_taxon = $db->get_taxon(-name => $f[1]);
    if (!defined $include_taxon) {
        print "    include taxon $f[1] not recognized\n";
        $anyWarnings++;
    }
    $groupings{$f[0]}{'include'} = $f[1];
    if ($f[2] ne "none") {
        my @excludeTaxa = split /\,/, $f[2];
        foreach my $excludeTaxon (@excludeTaxa) {
            ## check excludeTaxon exists in taxonomy database
            my $exclude_taxon = $db->get_taxon(-name => $excludeTaxon);
            if (!defined $exclude_taxon) {
                print "    exclude taxon $excludeTaxon not recognized\n";
                $anyWarnings++;
            } 
        }
        $groupings{$f[0]}{'exclude'} = \@excludeTaxa;
    }
}
close GROUPS;
if ($anyWarnings>0) {
    die "\n\nTerminating - there were unrecognized taxonomic terms in the groupings file\n\n";
}


## go through each file, work out common name
print "Going through each file\n";
foreach my $file (@ARGV) {
    my $species = $file; $species =~ s/\.txt$//; $species =~ s/^hits_//;
    $species =~ s/_/ /g;
    #if ($species ne "Ailuropoda melanoleuca") {next;} # xxxtemp
    #print "\n#### species $species\n";
    my $taxon = $db->get_taxon(-name => $species);
    
    my $lineageString = $tree_functions->get_lineage_string($taxon);
    my @l = split /\;/, $lineageString;
    
    ## go through each possible group and see if this lineageString is compatible
    my @goodGroups;
    foreach my $group (sort keys %groupings) {
        
        #### first check whether the INCLUDE group is present
        my $matched = 0; 
        my $excluded = 0;
        foreach my $l (@l) {
            if ($l eq $groupings{$group}{'include'}) {
                $matched = 1;
            }
        }
        #### if the INCLUDE group was matched, then we check for absence of the EXCLUDE matches
        if ($matched == 1) {
            if (defined $groupings{$group}{'exclude'}) {
                my $excludeRef = $groupings{$group}{'exclude'};
                my @exclude = @$excludeRef;
                foreach my $e (@exclude) {
                    foreach my $l (@l) {
                        if ($l eq $e) { $excluded = 1; }
                    }
                }
            }
            ## if we didn't see any of the excluded groups, it's good!
            if ($excluded == 0) { push @goodGroups, $group; }
        }
    }
    ###check we didn't match >1 group, and take action according to what group it is in (if any)
    my $finalGroup;
    if (@goodGroups == 0) { $finalGroup = "other"; } 
    if (@goodGroups > 1) {
        die "\n\nTerminating - species $species was in more than one group. goodGroups @goodGroups\n\nlineageString: $lineageString\n\n";
    } 
    if (@goodGroups == 1) {$finalGroup = $goodGroups[0];}
    print "species $species is in group $finalGroup\n";
    ## do something!
    if (!-e $finalGroup) { mkdir $finalGroup; }
    system("mv $file $finalGroup");

}
