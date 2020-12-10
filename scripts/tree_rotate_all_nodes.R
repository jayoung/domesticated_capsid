rotate_all_nodes <- function (tree, plotTree=FALSE, ...) {
    #changes the way a tree is displayed (reorders taxa from bottom to top)
    newtree<-tree
    num_nodes <- Nnode(tree)
    num_tips <- Ntip(tree)
    for (i in (num_tips+1):(num_tips+num_nodes)) {
        #cat("i",i,"\n")
        #check whether it's a polytomy
        num_descendents <- length(which(newtree$edge[,1]==i))
        #if not a polytomy
        if (num_descendents == 2) {newtree <- rotate(newtree,i)}
        else { #for polytomies
            cat ("Found a polytomy. num_descendents",num_descendents,"\n")
            for (j in 1:num_descendents) {
                other <- num_descendents + 1 - j
                cat (    "i",i,"j",j,"other",other,"\n")
                if (j >= other) {
                    cat ("Breaking the loop\n")
                    break
                }
                cat ("rotating j",j,"with other",other,"\n")
                newtree <- rotate(newtree,i,polytom=c(j,other))
            }
        }
    }
    if (plotTree) { plot.phylo(newtree, ...) }
    return(newtree)
}

### this version does not check for polytomies
rotate_all_nodes_2 <- function (tree, plotTree=FALSE, ...) {
    #changes the way a tree is displayed (reorders taxa from bottom to top)
    newtree<-tree
    num_nodes <- Nnode(tree)
    num_tips <- Ntip(tree)
    for (i in (num_tips+1):(num_tips+num_nodes)) {
        newtree <- rotate(newtree,i)
    }
    if (plotTree) { plot.phylo(newtree, ...) }
    return(newtree)
}

### this version checks for polytomies but does not rotate them
rotate_all_nodes_3 <- function (tree, plotTree=FALSE, ...) {
    #changes the way a tree is displayed (reorders taxa from bottom to top)
    newtree<-tree
    num_nodes <- Nnode(tree)
    num_tips <- Ntip(tree)
    for (i in (num_tips+1):(num_tips+num_nodes)) {
        #check whether it's a polytomy
        num_descendents <- length(which(newtree$edge[,1]==i))
        #if not a polytomy
        if (num_descendents == 2) {newtree <- rotate(newtree,i)}
    }
    if (plotTree) { plot.phylo(newtree, ...) }
    return(newtree)
}
