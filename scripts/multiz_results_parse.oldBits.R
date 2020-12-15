
## read input tables - old
inputs <- lapply(inputFiles, function(x) { 
    hasHeader <- c(FALSE,TRUE)
    names(hasHeader) <- c("mafParse","pseudReport")
    out <- lapply(names(x), function(y) {
        read.delim(paste(multizDir, x[[y]], sep="/"), header=hasHeader[y]) 
    })
    names(out) <- names(x)
    return(out)
})


## convert those to a table of what's going on in each species in the tree
multiz_results_old <- lapply(names(inputs), readMultizResults_old)
names(multiz_results_old) <- names(inputs)


readMultizResults_old <- function(analysisName, refAssembly="hg38", 
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
