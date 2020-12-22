setwd("/Users/jayoung/Desktop/mac_workStuff/domesticated_capsid")
library(ape)
source("scripts/tree_rotate_all_nodes.R")
source("scripts/multiz_results_functions.R")

######### read in species info, tree, etc

species_dat <- read.delim("/Users/jayoung/Desktop/mac_workStuff/lib/hg38.README.tabdelim.txt")
# make sure I know which are mammals
species_dat[,"isMammal"] <- ! species_dat[,"Species.set"] %in% c("Aves","Fish","Sarcopterygii")
# table(species_dat[,"isMammal"])
#  FALSE  TRUE 
#     38    62 

#### species tree
tree <- read.tree("/Users/jayoung/Desktop/mac_workStuff/lib/hg38.100way.nh")
## just mammals
tree_mammals <- keep.tip(tree, species_dat[which(species_dat[,"isMammal"]),"Assembly.name"])
## just placental mammals
tree_placMammals <- drop.tip(tree_mammals, c("monDom5","sarHar1","macEug2","ornAna1"))
## a version with common species names
tree_placMammals_commonNames <- tree_placMammals
tree_placMammals_commonNames$tip.label <- species_dat[ match(tree_placMammals_commonNames$tip.label, species_dat[,"Assembly.name"]), "Common.name" ]   
tree_placMammals_commonNames <- rotate_all_nodes(tree_placMammals_commonNames)


##### read in info from parsed multiz alignments

multizDir <- "multiz_alignments"

mafLogFiles <- list.files(multizDir, pattern=".maf2fa.log.txt$")
pseudReportFiles <- list.files(multizDir, pattern=".pseudReportEachGene.txt$")

analysisNames <- gsub(".hg38.100way_placMamm.maf2fa.log.txt","",mafLogFiles)

multiz_results <- lapply(analysisNames, readMultizResults)
names(multiz_results) <- analysisNames

summaryTableByGene <- t(sapply(multiz_results, function(x) { table(x[,"status"])}))
#                           Reference Intact Truncated Pseud Absent
# ARC_NM_015193                     1     43         2     6      6
# CCDC8_NM_032040                   1      9         1    45      2
# LDOC1_NM_012317                   1     52         1     3      1
# etc...

## check species are in the same order
table( sapply(2:length(multiz_results), function(i) {
    identical(rownames(multiz_results[[1]]), rownames(multiz_results[[i]]))
} ) )
# TRUE 
#   45 

tableBySpecies <- sapply(multiz_results, function(x) { x[,"status"] })
rownames(tableBySpecies) <- rownames(multiz_results[[1]])

summaryBySpecies <- t(apply(tableBySpecies, 1, function(x) {
    table(factor(x, levels=c("Reference", "Intact", "Truncated", "Pseud", "Absent")))
}))


pdf(height=7,width=11, file="plots/geneStatus_plots_multizMafs.pdf")
plotStatusManyGenes(analysisNames)
dev.off()





##### read in info from parsed liftOver alignments

liftOver_pseudReportFiles <- list.files("liftOver_alignments/liftOverFiles_justCDS/alignments", pattern=".pseudReportEachGene.txt$")
names(liftOver_pseudReportFiles) <- gsub(".cdsSeqs_aln1_NT.sorted.pseudReportEachGene.txt","", liftOver_pseudReportFiles)

liftOverLogFiles <- list.files("liftOver_alignments/liftOverFiles_justCDS/liftOverFiles", pattern=".getCDSlog.txt", full.names=TRUE)
names(liftOverLogFiles) <- sapply(strsplit(liftOverLogFiles, "\\."), "[[", 3)

liftOverLogs <- lapply(liftOverLogFiles, scan, what="character", sep="\n", quiet=TRUE)

## there are some species for which liftOver didn't work at all - figure out which ones those are, and make a new tree to use for plots etc
liftOverImpossible <- names(liftOverLogs) [which( sapply(liftOverLogs, function(x) {
    sum(grepl("Assembly not supported by ucscApiClient", x))>0
}) )]

## get a tree with only the species where liftOver is possible:
tree_placMammals_liftOverSubset <- drop.tip(tree_placMammals, liftOverImpossible)
## a version with common species names
tree_placMammals_liftOverSubset_commonNames <- tree_placMammals_liftOverSubset
tree_placMammals_liftOverSubset_commonNames$tip.label <- species_dat[ match(tree_placMammals_liftOverSubset_commonNames$tip.label, species_dat[,"Assembly.name"]), "Common.name" ]   
tree_placMammals_liftOverSubset_commonNames <- rotate_all_nodes(tree_placMammals_liftOverSubset_commonNames)


liftOver_results <- lapply(liftOver_pseudReportFiles, 
                           readLiftOverResults, 
                           dir="liftOver_alignments/liftOverFiles_justCDS/alignments")

summaryTableByGene_liftOver <- t(sapply(liftOver_results, function(x) { table(x[,"status"])}))
#                           Reference Intact Truncated Pseud Absent
# ARC_NM_015193                     1     37         0     5     15
# CCDC8_NM_032040                   1     37         6     7      7
# LDOC1_NM_012317                   1     47         0     3      7
# etc...

## check species are in the same order
table( sapply(2:length(liftOver_results), function(i) {
    identical(rownames(liftOver_results[[1]]), rownames(liftOver_results[[i]]))
} ) )
# TRUE 
#   45 

tableBySpecies_liftOver <- sapply(liftOver_results, function(x) { x[,"status"] })
rownames(tableBySpecies_liftOver) <- rownames(liftOver_results[[1]])

summaryBySpecies_liftOver <- t(apply(tableBySpecies_liftOver, 1, function(x) {
    table(factor(x, levels=c("Reference", "Intact", "Truncated", "Pseud", "Absent")))
}))


pdf(height=7,width=11, file="plots/geneStatus_plots_liftOverAlns_justCDS.pdf")
plotStatusManyGenes(analysisNames, 
                    tree=tree_placMammals_liftOverSubset_commonNames,
                    resultsTables=liftOver_results,
                    myTitle="Gene status in placental mammals\n(liftOver seqs, just CDS, MACSE aligned)")
dev.off()
