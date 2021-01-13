setwd("/Users/jayoung/Desktop/mac_workStuff/domesticated_capsid")
library(ape)
library(openxlsx)
source("scripts/tree_rotate_all_nodes.R")
source("scripts/parseAlignmentReports_functions.R")


######### read in species info
species_dat <- read.delim("miscellaneous/hg38.README.tabdelim.addSpecies.txt")
# make sure I know which are mammals
species_dat[,"isMammal"] <- ! species_dat[,"Species.set"] %in% c("Aves","Fish","Sarcopterygii")
# table(species_dat[,"isMammal"])
#  FALSE  TRUE 
#     39    63
species_dat[,"Assembly.name.short"] <- gsub("\\d+$","",species_dat[,"Assembly.name"], perl=TRUE)

######### read in species tree

### I copied this from /Volumes/malik_h/grp/public_databases/UCSC/human_Dec2013/conservation_tracks/hg38.100way.scientificNames.nh 
# and manually added species that were missing (i.e. species I included from blast that were not in the 100-way set)

tree_latinNames <- read.tree("miscellaneous/hg38.100way.scientificNames.addSpecies.nh")

## get tree of just mammals
speciesToKeep <- gsub(" ", "_",species_dat[which(species_dat[,"isMammal"]),"Latin.name"])
# special case, where there was a mismatch between UCSC tree file and species info file"
speciesToKeep <- gsub("hamadryas","anubis",speciesToKeep)
tree_latinNames_mammals <- keep.tip(tree_latinNames, speciesToKeep)

## a version with common species names - not sure if I need this
tree_commonNames <- tree_latinNames
tree_commonNames$tip.label <- species_dat[ match(gsub("_", " ",tree_commonNames$tip.label),
                                        species_dat[,"Latin.name"]), "Common.name" ]   
tree_commonNames <- rotate_all_nodes(tree_commonNames)


##### read in info from parsed blast-based alignments

blastAlnsDirs <- list.files(list.files("alignments", full.names=TRUE), pattern="masterAlignments", full.names=TRUE)

pseudReportFiles <- list.files(blastAlnsDirs, pattern=".pseudReportEachGene.txt$", full.names=TRUE)

analysisNames <- sapply( strsplit(pseudReportFiles, "/"), function(x) {x[length(x)]})
analysisNames <- sapply( strsplit(analysisNames, "\\."), function(x) { 
    if (grepl("^ARC_", x)[1]) {
        y <- "ARC_ARC"
    } else {
        y <- grep("^group_", x, value=TRUE)[1] 
        y <- gsub("^group_", "", y)
        y <- gsub("_aln\\d+_NT", "", y)
    }
    return(y)
})
    
# xxx need to deal with situations where there is >1 sequence for a species in an alignment e.g pseudReports[["PNMA_PNMA6EF"]]


pseudReports <- lapply(pseudReportFiles, readBlastBasedPseudReport)
names(pseudReports) <- analysisNames

# simple list of species in any analysis
speciesAnalyzed <- unique(  unlist( sapply(pseudReports, "[[", "Species"))  )
speciesAnalyzed_tree <- rotate_all_nodes(keep.tip(tree_latinNames, gsub(" ","_",speciesAnalyzed)))

# counts by species 
pseudReportsCounts <- lapply(pseudReports, function(x) {
    x[,"Species"] <- factor(x[,"Species"], levels=gsub("_"," ",speciesAnalyzed_tree$tip.label))
    y <- table(x[,"Species"], x[,"Pseud"] )
    ## for species that weren't in the pseudReport, replace 0 with NA
    speciesNotAnalyzed <- setdiff(rownames(y), unique(x[,"Species"]))
    y[speciesNotAnalyzed, ] <- NA
    return(y)
})

## save in an Excel file
saveReports(pseudReportsCounts, "alignments/reports/pseudReports_blastBasedAlignments_2021_Jan13.xlsx")

pdf(height=7,width=11, file="plots/pseudReports_blastBasedAlignments_Summary_2021_Jan13.pdf")
plot.phylo(speciesAnalyzed_tree, 
           cex=1, font=1, 
           align.tip.label=TRUE,
           x.lim=3, 
           no.margin = TRUE) 
dev.off()


#xxxx old code:
pdf(height=7,width=11, file="plots/geneStatus_plots_multizMafs.pdf")
plotStatusManyGenes(analysisNames)
dev.off()





pdf(height=7,width=11, file="plots/geneStatus_plots_liftOverAlns_justCDS.pdf")
plotStatusManyGenes(analysisNames, 
                    tree=tree_placMammals_liftOverSubset_commonNames,
                    resultsTables=liftOver_results,
                    myTitle="Gene status in placental mammals\n(liftOver seqs, just CDS, MACSE aligned)")
dev.off()
