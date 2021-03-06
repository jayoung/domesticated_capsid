setwd("/Users/jayoung/Desktop/mac_workStuff/domesticated_capsid")
library(ape)
source("scripts/tree_rotate_all_nodes.R")
source("scripts/parseAlignmentReports_functions.R")

######### read in species info, tree, etc

species_dat <- read.delim("/Users/jayoung/Desktop/mac_workStuff/lib/hg38.README.tabdelim.txt")
# make sure I know which are mammals
species_dat[,"isMammal"] <- ! species_dat[,"Species.set"] %in% c("Aves","Fish","Sarcopterygii")
# table(species_dat[,"isMammal"])
#  FALSE  TRUE 
#     38    62 

species_dat[,"Assembly.name.short"] <- gsub("\\d+$","",species_dat[,"Assembly.name"], perl=TRUE)

### get newer assembly names:
newerAssemblies <- scan("/Users/jayoung/Desktop/mac_workStuff/lib/hg38.README.placentalMammals.newerAssemblies.txt", what="character", sep="\n")
temp <- gsub("\\d+$","",newerAssemblies, perl=TRUE)
table(temp %in% species_dat[,"Assembly.name.short"])
species_dat[,"newerAssembly"] <- NA
species_dat[match(temp, species_dat[,"Assembly.name.short"]),"newerAssembly"] <- newerAssemblies
rm(temp)

#### species tree
tree <- read.tree("/Users/jayoung/Desktop/mac_workStuff/lib/hg38.100way.nh")
## just mammals
tree_mammals <- keep.tip(tree, species_dat[which(species_dat[,"isMammal"]),"Assembly.name"])

## just placental mammals
tree_placMammals <- drop.tip(tree_mammals, c("monDom5","sarHar1","macEug2","ornAna1"))

tree_placMammals_shortAssemblyNames <- tree_placMammals
tree_placMammals_shortAssemblyNames$tip.label <- gsub("\\d+$","",tree_placMammals_shortAssemblyNames$tip.label, perl=TRUE)

tree_placMammals_newAssemblyNames <- tree_placMammals_shortAssemblyNames
tree_placMammals_newAssemblyNames$tip.label <- species_dat[match(tree_placMammals_shortAssemblyNames$tip.label, species_dat[,"Assembly.name.short"]),"newerAssembly"]


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





##### read in info from parsed liftOver alignments, just CDS, same assemblies as for 100-way mafs

liftOver_pseudReportFiles <- list.files("liftOver_alignments/liftOverFiles_justCDS/alignments", pattern=".pseudReportEachGene.txt$")
names(liftOver_pseudReportFiles) <- gsub(".cdsSeqs_aln1_NT.sorted.pseudReportEachGene.txt","", liftOver_pseudReportFiles)

liftOverLogFiles <- list.files("liftOver_alignments/liftOverFiles_justCDS/liftOverFiles", pattern=".getCDSlog.txt", full.names=TRUE)
names(liftOverLogFiles) <- sapply(strsplit(liftOverLogFiles, "\\."), "[[", 5)

liftOverLogs <- lapply(liftOverLogFiles, scan, what="character", sep="\n", quiet=TRUE)

## there are some species for which liftOver didn't work at all - figure out which ones those are, and make a new tree to use for plots etc
liftOverImpossible <- names(liftOverLogs) [which( sapply(liftOverLogs, function(x) {
    sum(grepl("Assembly not supported by ucscApiClient", x))>0
}) )]


liftOverSeqsWithNs <- lapply( liftOverLogs, function(x) { 
    y <- grep("of bases are Ns", x, value=TRUE)
    if (length(y)==0) {return(y)}
    y <- gsub("    WARNING - ignoring seq ","",y)
    y <- sapply(strsplit(y, ","), "[[", 1)
    assembly <- strsplit(y[1], "_")[[1]]
    assembly <- assembly[length(assembly)]
    y <- gsub( paste("_",assembly,sep=""),"",y)
    y <- sapply(strsplit(y, "\\."), "[[", 1)
    return(y)
})


## get a tree with only the species where liftOver is possible:
tree_placMammals_liftOverSubset <- drop.tip(tree_placMammals, liftOverImpossible)
## a version with common species names
tree_placMammals_liftOverSubset_commonNames <- tree_placMammals_liftOverSubset
tree_placMammals_liftOverSubset_commonNames$tip.label <- species_dat[ match(tree_placMammals_liftOverSubset_commonNames$tip.label, species_dat[,"Assembly.name"]), "Common.name" ]   
tree_placMammals_liftOverSubset_commonNames <- rotate_all_nodes(tree_placMammals_liftOverSubset_commonNames)


liftOver_results <- lapply(liftOver_pseudReportFiles, 
                           readLiftOverResults, 
                           dir="liftOver_alignments/liftOverFiles_justCDS/alignments")

## note which seqs were truncated by Ns so much that I didn't try to align them (from liftOverSeqsWithNs)
liftOver_results <- lapply( analysisNames, function(x)  {
    geneName <- strsplit(x, "_")[[1]][1]
    transcriptName <- gsub( paste(geneName,"_",sep=""), "", x)
    
    speciesWithNs <- names(liftOverSeqsWithNs)[grep(transcriptName, liftOverSeqsWithNs)]
    #return(speciesWithNs)
    results <- liftOver_results[[x]]
    results[match(speciesWithNs, rownames(results)),"status"] <- "Truncated"
    return(results)
    
} )
names(liftOver_results) <- analysisNames


summaryTableByGene_liftOver <- t(sapply(liftOver_results, function(x) { table(x[,"status"])}))
head(summaryTableByGene_liftOver)
#                           Reference Intact Truncated Pseud Absent
# ARC_NM_015193                     1     37         2     5     13
# CCDC8_NM_032040                   1     37        10     7      3
# LDOC1_NM_012317                   1     47         1     3      6

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




##### read in info from parsed liftOver alignments, just CDS, newer assemblies

liftOver_pseudReportFiles <- list.files("liftOver_alignments/liftOverFiles_justCDS_newAssemblies/alignments", pattern=".pseudReportEachGene.txt$")
names(liftOver_pseudReportFiles) <- gsub(".cdsSeqs_aln1_NT.sorted.pseudReportEachGene.txt","", liftOver_pseudReportFiles)

liftOverLogFiles <- list.files("liftOver_alignments/liftOverFiles_justCDS_newAssemblies/liftOverFiles", pattern=".getCDSlog.txt", full.names=TRUE)
names(liftOverLogFiles) <- sapply(strsplit(liftOverLogFiles, "\\."), "[[", 5)

liftOverLogs <- lapply(liftOverLogFiles, scan, what="character", sep="\n", quiet=TRUE)

## there are some species for which liftOver didn't work at all - figure out which ones those are, and make a new tree to use for plots etc
liftOverImpossible <- names(liftOverLogs) [which( sapply(liftOverLogs, function(x) {
    sum(grepl("Assembly not supported by ucscApiClient", x))>0
}) )]
# there were none

liftOverSeqsWithNs <- lapply( liftOverLogs, function(x) { 
    y <- grep("of bases are Ns", x, value=TRUE)
    if (length(y)==0) {return(y)}
    y <- gsub("    WARNING - ignoring seq ","",y)
    y <- sapply(strsplit(y, ","), "[[", 1)
    assembly <- strsplit(y[1], "_")[[1]]
    assembly <- assembly[length(assembly)]
    y <- gsub( paste("_",assembly,sep=""),"",y)
    y <- sapply(strsplit(y, "\\."), "[[", 1)
    return(y)
})



# xxx get a version of tree_placMammals that has short names

liftOver_results <- lapply(liftOver_pseudReportFiles, 
                    readLiftOverResults, 
                    assemblyColNameInfo="newerAssembly",
                    refAssembly="hg38",
                    dir="liftOver_alignments/liftOverFiles_justCDS_newAssemblies/alignments",
                    tree=tree_placMammals_newAssemblyNames)

## note which seqs were truncated by Ns so much that I didn't try to align them (from liftOverSeqsWithNs)
liftOver_results <- lapply( analysisNames, function(x)  {
    geneName <- strsplit(x, "_")[[1]][1]
    transcriptName <- gsub( paste(geneName,"_",sep=""), "", x)
    
    speciesWithNs <- names(liftOverSeqsWithNs)[grep(transcriptName, liftOverSeqsWithNs)]
    #return(speciesWithNs)
    results <- liftOver_results[[x]]
    results[match(speciesWithNs, rownames(results)),"status"] <- "Truncated"
    return(results)
    
} )
names(liftOver_results) <- analysisNames


summaryTableByGene_liftOver <- t(sapply(liftOver_results, function(x) { table(x[,"status"])}))
head(summaryTableByGene_liftOver)
#                           Reference Intact Truncated Pseud Absent
# ARC_NM_015193                     1     37         2     5     13
# CCDC8_NM_032040                   1     37        10     7      3
# LDOC1_NM_012317                   1     47         1     3      6

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
