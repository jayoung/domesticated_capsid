readMultizResults <- function(analysisName, 
                              resultsDir="multiz_alignments",
                              refAssembly="hg38", 
                              tree=tree_placMammals, 
                              info=species_dat,
                              mafParseFileSuffix=".hg38.100way_placMamm.maf2fa.log.txt",
                              pseudReportFileSuffix=".hg38.100way_placMamm.fa.nogaps_aln1_NT.sorted.pseudReportEachGene.txt") {
    
    ## start the table using the tree and info table
    results <- data.frame(assemblyName=tree$tip.label, row.names=tree$tip.label)
    results[,"commonName"] <- info[ match(results[,"assemblyName"], info[,"Assembly.name"]), "Common.name" ]
    results[,"status"] <- factor(NA, 
                                 levels=c("Reference", "Intact", "Truncated", "Pseud", "Absent"))
    results[which(results[,"assemblyName"]==refAssembly),"status"] <- "Reference"

    ## read in the mafParsing log and the pseudogene report
    mafParseFile <- paste(resultsDir, "/", analysisName, mafParseFileSuffix, sep="")
    pseudReportFile <- paste(resultsDir, "/", analysisName, pseudReportFileSuffix, sep="")
    if(!file.exists(mafParseFile)) {
        stop("ERROR - cannot find mafParseFile",mafParseFile," - check mafParseFileSuffix option\n")
    }
    if(!file.exists(pseudReportFile)) {
        stop("ERROR - cannot find pseudReportFile",pseudReportFile," - check pseudReportFileSuffix option\n")
    }
    mafParse <- read.delim(mafParseFile, header=FALSE) 
    pseudReport <- read.delim(pseudReportFile, header=TRUE) 

    ## add the seqs that were missing from maf using x[["mafParse"]]
    missingFromMaf <- mafParse
    missingFromMaf <- missingFromMaf[which(missingFromMaf[,1] != refAssembly),]
    missingFromMaf <- missingFromMaf[,1]
    results[which(results[,"assemblyName"] %in% missingFromMaf),"status"] <- "Absent"
    
    ## add intact/pseud status using x[["pseudReport"]]
    pseuds <- pseudReport
    pseuds <- pseuds[which(pseuds[,"Pseud"] != "Reference"),]
    results[ match(pseuds[,"Seq"], results[,"assemblyName"]), "status"] <- pseuds[,"Pseud"]
    return(results)
} 




readLiftOverResults <- function(file=liftOver_pseudReportFile, 
                                dir="liftOver_alignments/alignments",
                                refAssembly="hg38", 
                                tree=tree_placMammals_liftOverSubset, 
                                info=species_dat,
                                pseudReportFileSuffix=".cdsSeqs_aln1_NT.sorted.pseudReportEachGene.txt") {
    
    ## start the table using the tree and info table
    results <- data.frame(assemblyName=tree$tip.label, row.names=tree$tip.label)
    results[,"commonName"] <- info[ match(results[,"assemblyName"], info[,"Assembly.name"]), "Common.name" ]
    results[,"status"] <- factor(NA, 
                                 levels=c("Reference", "Intact", "Truncated", "Pseud", "Absent"))
    results[which(results[,"assemblyName"]==refAssembly),"status"] <- "Reference"
    
    ## read in the pseudogene report
    if(!file.exists(paste(dir,file,sep="/"))) {
        stop("ERROR - cannot find pseudReportFile",paste(dir,file,sep="/")," - check pseudReportFileSuffix option\n")
    }
    pseudReport <- read.delim(paste(dir,file,sep="/"), header=TRUE) 
    pseudReport[,"species"] <- sapply( strsplit(pseudReport[,"Seq"],"_"), function(x) {x[length(x)]})
    
    ## add the seqs that were missing from the pseudReport (they failed to liftOver)
    results[which(! results[,"assemblyName"] %in% pseudReport[,"species"]),"status"] <- "Absent"
    
    ## add intact/pseud status using x[["pseudReport"]]
    pseuds <- pseudReport
    pseuds <- pseuds[which(pseuds[,"Pseud"] != "Reference"),]
    results[ match(pseuds[,"species"], results[,"assemblyName"]), "status"] <- pseuds[,"Pseud"]
    return(results)
} 


#### define how I want to plot each type of gene
geneStatusSymbols <- list()
geneStatusSymbols[["Reference"]][["pch"]] <- 15
geneStatusSymbols[["Intact"]][["pch"]] <- 15
geneStatusSymbols[["Pseud"]][["pch"]] <- 7
geneStatusSymbols[["Truncated"]][["pch"]] <- 22
geneStatusSymbols[["Truncated"]][["bg"]] <- "gray"
geneStatusSymbols[["Absent"]][["pch"]] <- 0

#### add intact/pseud boxes to an already-plotted tree
plotGeneStatusSymbols <- function(analysisName, Xoffset=0,
                                  tree=tree_placMammals_commonNames,
                                  resultsTables=multiz_results,
                                  key=geneStatusSymbols, 
                                  myTitle=NULL) {
    ### get X and Y positions of terminal nodes 
    numSpec <- length(tree$tip.label)
    alignedXpos <- max(node.depth.edgelength(tree))
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    yPos <- lastPP$yy[1:numSpec] ## y positions of terminal nodes
    
    ##### get pch and bg for each gene (plot characters)
    res <- resultsTables[[analysisName]]
    res[,"pch"] <- NA
    res[,"bg"] <- "black"
    for (geneClass in names(key)) {
        res[which(res[,"status"]==geneClass),"pch"] <- key[[geneClass]][["pch"]]
        if("bg" %in% names(key[[geneClass]])) {
            res[which(res[,"status"]==geneClass),"bg"] <- key[[geneClass]][["bg"]]
        }
    }
    myPch <- res[match(tree$tip.label, res[,"commonName"]), "pch"]
    myBg <- res[match(tree$tip.label, res[,"commonName"]), "bg"]
    points(x=rep(alignedXpos+Xoffset, numSpec), 
           y=yPos,
           pch=myPch, cex = 0.75, lwd=0.5, bg=myBg)
    if(!is.null(myTitle)) {
        text( x=alignedXpos+Xoffset, 
              y=max(yPos) + 1, adj=0,
              myTitle, cex=0.5, srt=90) 
    }
}

plotGeneStatusLegend <- function(pos, key=geneStatusSymbols, ...) {
    keyText <- names(key)
    keyPch <- sapply(key, "[[", "pch")
    keyBg <- lapply(key, "[[", "bg")
    keyBg[which(sapply(keyBg, is.null))] <- "black"
    keyBg <- unlist(keyBg)
    legend(x=pos, legend=keyText, pch=keyPch, pt.bg=keyBg, ...)
}


plotStatusManyGenes <- function(analyses, 
                                tree=tree_placMammals_commonNames, 
                                resultsTables=multiz_results,
                                myTitle="Gene status in placental mammals\n(100-way mafs, MACSE realigned)",
                                xSeparation=0.02,
                                xAmountToAddForSpeciesNames=0.51,
                                yAmountToAddForIsoformNames=25) {
    ## get num genes so I can put a space between each gene, add some more
    numGenes <- length(unique(sapply(strsplit(analyses, "_"), "[[", 1)))
    tempLabelOffset <- xSeparation/2 + (length(analyses)+numGenes-1) * xSeparation
    cat("num genes",numGenes,"\n")
    plot.phylo(tree, 
               cex=0.5, font=1, 
               label.offset = tempLabelOffset,
               align.tip.label=TRUE,
               x.lim=tempLabelOffset+xAmountToAddForSpeciesNames, 
               y.lim=length(tree$tip.label) + 
                            yAmountToAddForIsoformNames, 
               no.margin = TRUE) 
    currentOffset <- 0.01
    prevGene <- strsplit(analyses[1], "_")[[1]][1]
    for(thisAnalysis in analyses) {
        geneName <- strsplit(thisAnalysis, "_")[[1]][1]
        if (geneName != prevGene) { 
            currentOffset <- currentOffset + xSeparation 
        }
        prevGene <- geneName
        plotGeneStatusSymbols(thisAnalysis, 
                              tree=tree,
                              Xoffset=currentOffset, 
                              resultsTables=resultsTables,
                              myTitle=thisAnalysis)
        currentOffset <- currentOffset + xSeparation
    }
    plotGeneStatusLegend("topleft", cex=1)
    title(myTitle, line=-2, cex.main=1)
}

