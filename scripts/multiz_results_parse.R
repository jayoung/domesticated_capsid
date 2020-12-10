library(ape)
source("/Users/jayoung/Desktop/mac_workStuff/mac_domesticated_gags_BlakeWiedenheft/lib/tree_rotate_all_nodes.R")
source("multiz_results_functions.R")

#### read in species info
species_dat <- read.delim("/Users/jayoung/Desktop/mac_workStuff/mac_domesticated_gags_BlakeWiedenheft/lib/hg38.README.tabdelim.txt")

# make sure I know which are mammals
species_dat[,"isMammal"] <- ! species_dat[,"Species.set"] %in% c("Aves","Fish","Sarcopterygii")
# table(species_dat[,"isMammal"])
#  FALSE  TRUE 
#     38    62 

##### read in species tree
tree <- read.tree("/Users/jayoung/Desktop/mac_workStuff/mac_domesticated_gags_BlakeWiedenheft/lib/hg38.100way.nh")

## just mammals
tree_mammals <- keep.tip(tree, species_dat[which(species_dat[,"isMammal"]),"Assembly.name"])

## just placental mammals
tree_placMammals <- drop.tip(tree_mammals, c("monDom5","sarHar1","macEug2","ornAna1"))

## a version with common species names
tree_placMammals_commonNames <- tree_placMammals
tree_placMammals_commonNames$tip.label <- species_dat[ match(tree_placMammals_commonNames$tip.label, species_dat[,"Assembly.name"]), "Common.name" ]   

tree_placMammals_commonNames <- rotate_all_nodes(tree_placMammals_commonNames)


##### read in some of the parsed multiz stuff

inputFiles <- list()
inputFiles[["PNMA3_NM_013364"]] <- list()
inputFiles[["PNMA3_NM_013364"]][["mafParse"]] <- "PNMA3_NM_013364.refGene.100way.frags.2.combined.log.txt"
inputFiles[["PNMA3_NM_013364"]][["pseudReport"]] <- "PNMA3_NM_013364.refGene.100way.frags.2.combined.fa.nogaps.MACSE.nt.sorted.pseudReportEachGene.txt"

inputFiles[["PNMA3_NM_001282535"]] <- list()
inputFiles[["PNMA3_NM_001282535"]][["mafParse"]] <- "PNMA3_NM_001282535.refGene.100way.frags.2.combined.log.txt"
inputFiles[["PNMA3_NM_001282535"]][["pseudReport"]] <- "PNMA3_NM_001282535.refGene.100way.frags.2.combined.fa.nogaps.MACSE.nt.sorted.pseudReportEachGene.txt"

## read input tables
inputs <- lapply(inputFiles, function(x) { 
    hasHeader <- c(FALSE,TRUE)
    names(hasHeader) <- c("mafParse","pseudReport")
    out <- lapply(names(x), function(y) {
        read.delim(x[[y]], header=hasHeader[y]) 
    })
    names(out) <- names(x)
    return(out)
})

## convert those to a table of what's going on in each species in the tree
multiz_results <- lapply(names(inputs), readMultizResults)
names(multiz_results) <- names(inputs)

t(sapply(multiz_results, function(x) { table(x[,"status"])}))
#                    Absent Intact Pseud Reference Truncated
# PNMA3_NM_013364        19     16    21         1         1
# PNMA3_NM_001282535     18     11    26         1         2



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



