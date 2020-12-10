
readMultizResults <- function(analysisName, refAssembly="hg38", 
                              tree=tree_placMammals, info=species_dat) {
    
    ## start the table using the tree and info table
    results <- data.frame(assemblyName=tree$tip.label, row.names=tree$tip.label)
    results[,"commonName"] <- info[ match(results[,"assemblyName"], info[,"Assembly.name"]), "Common.name" ]
    results[,"status"] <- NA
    results[which(results[,"assemblyName"]==refAssembly),"status"] <- "Reference"
    
    ## add the seqs that were missing from maf using x[["mafParse"]]
    missingFromMaf <- inputs[[analysisName]][["mafParse"]]
    missingFromMaf <- missingFromMaf[which(missingFromMaf[,1] != refAssembly),]
    missingFromMaf <- missingFromMaf[,1]
    results[which(results[,"assemblyName"] %in% missingFromMaf),"status"] <- "Absent"
    
    ## add intact/pseud status using x[["pseudReport"]]
    pseuds <- inputs[[analysisName]][["pseudReport"]]
    pseuds <- pseuds[which(pseuds[,"Pseud"] != "Reference"),]
    results[ match(pseuds[,"Seq"], results[,"assemblyName"]), "status"] <- pseuds[,"Pseud"]
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


