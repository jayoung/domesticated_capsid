setwd("/Users/jayoung/Desktop/mac_workStuff/domesticated_capsid")

library(ape)
source("scripts/tree_rotate_all_nodes.R")
source("scripts/multiz_results_functions.R")

#### read in species info
species_dat <- read.delim("/Users/jayoung/Desktop/mac_workStuff/lib/hg38.README.tabdelim.txt")

# make sure I know which are mammals
species_dat[,"isMammal"] <- ! species_dat[,"Species.set"] %in% c("Aves","Fish","Sarcopterygii")
# table(species_dat[,"isMammal"])
#  FALSE  TRUE 
#     38    62 

##### read in species tree
tree <- read.tree("/Users/jayoung/Desktop/mac_workStuff/lib/hg38.100way.nh")

## just mammals
tree_mammals <- keep.tip(tree, species_dat[which(species_dat[,"isMammal"]),"Assembly.name"])

## just placental mammals
tree_placMammals <- drop.tip(tree_mammals, c("monDom5","sarHar1","macEug2","ornAna1"))

## a version with common species names
tree_placMammals_commonNames <- tree_placMammals
tree_placMammals_commonNames$tip.label <- species_dat[ match(tree_placMammals_commonNames$tip.label, species_dat[,"Assembly.name"]), "Common.name" ]   

tree_placMammals_commonNames <- rotate_all_nodes(tree_placMammals_commonNames)


##### read in some of the parsed multiz stuff

multizDir <- "multiz_alignments"

## old:
inputFiles <- list()
inputFiles[["PNMA3_NM_013364"]] <- list()
inputFiles[["PNMA3_NM_013364"]][["mafParse"]] <- "PNMA3_NM_013364.hg38.100way_placMamm.maf2fa.log.txt"
inputFiles[["PNMA3_NM_013364"]][["pseudReport"]] <- "PNMA3_NM_013364.hg38.100way_placMamm.fa.nogaps_aln1_NT.sorted.pseudReportEachGene.txt"

inputFiles[["PNMA3_NM_001282535"]] <- list()
inputFiles[["PNMA3_NM_001282535"]][["mafParse"]] <- "PNMA3_NM_001282535.hg38.100way_placMamm.maf2fa.log.txt"
inputFiles[["PNMA3_NM_001282535"]][["pseudReport"]] <- "PNMA3_NM_001282535.hg38.100way_placMamm.fa.nogaps_aln1_NT.sorted.pseudReportEachGene.txt"



## new:

mafLogFiles <- list.files(multizDir, pattern=".maf2fa.log.txt$")
pseudReportFiles <- list.files(multizDir, pattern=".pseudReportEachGene.txt$")

source("scripts/multiz_results_functions.R")

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


## test plot, just two isoforms
plot.phylo(tree_placMammals_commonNames, 
                     cex=0.5, font=1, label.offset = 0.05,
                     align.tip.label=TRUE,
                     x.lim=0.55, 
                     y.lim=length(tree_placMammals_commonNames$tip.label)+17, 
                     no.margin = TRUE) 
plotGeneStatusSymbols("PNMA3_NM_013364", Xoffset=0.01, myTitle="PNMA3_NM_013364")
plotGeneStatusSymbols("PNMA3_NM_001282535", Xoffset=0.03, myTitle="PNMA3_NM_001282535")
plotGeneStatusLegend("topleft", cex=1)
title("Gene status in placental mammals\n(100-way mafs, MACSE realigned)", line=-2, cex.main=1, adj=1)


plotStatusManyGenes <- function(analyses, 
                                xSeparation=0.02,
                                xAmountToAddForSpeciesNames=0.51,
                                yAmountToAddForIsoformNames=25) {
    ## get num genes so I can put a space between each gene, add some more
    numGenes <- length(unique(sapply(strsplit(analyses, "_"), "[[", 1)))
    tempLabelOffset <- xSeparation/2 + (length(analyses)+numGenes-1) * xSeparation
    cat("num genes",numGenes,"\n")
    plot.phylo(tree_placMammals_commonNames, 
               cex=0.5, font=1, 
               label.offset = tempLabelOffset,
               align.tip.label=TRUE,
               x.lim=tempLabelOffset+xAmountToAddForSpeciesNames, 
               y.lim=length(tree_placMammals_commonNames$tip.label)+yAmountToAddForIsoformNames, 
               no.margin = TRUE) 
    currentOffset <- 0.01
    prevGene <- strsplit(analyses[1], "_")[[1]][1]
    for(thisAnalysis in analyses) {
        geneName <- strsplit(thisAnalysis, "_")[[1]][1]
        if (geneName != prevGene) { 
            currentOffset <- currentOffset + xSeparation 
        }
        prevGene <- geneName
        plotGeneStatusSymbols(thisAnalysis, Xoffset=currentOffset, myTitle=thisAnalysis)
        currentOffset <- currentOffset + xSeparation
    }
    plotGeneStatusLegend("topleft", cex=1)
    title("Gene status in placental mammals\n(100-way mafs, MACSE realigned)", 
          line=-2, cex.main=1)
}
pdf(height=7,width=11, file="plots/multiz_status_plots.pdf")
plotStatusManyGenes(analysisNames)
dev.off()
