.checkphyloarg <- function(phyl){
  if(inherits(phyl,"phylo4d"))
    phyl <- extractTree(phyl)
  if(inherits(phyl,"phylo4"))    
    phyl.phylo <- as(phyl,"phylo")
  if(inherits(phyl,"phylo")){
    phyl.phylo <- phyl
    phyl <- as(phyl.phylo,"phylo4")
  }
  if(inherits(phyl, "hclust")){
    phyl.phylo <- as.phylo(phyl)
    phyl <- as(phyl.phylo,"phylo4")
}
  if(!exists("phyl.phylo"))
    stop("unconvenient phyl: must be a hclust, phylo, phylog, phylo4 or phylo4d object")

  if(!isRooted(phyl)){
    phyl.phylo$root.edge <- 0
    phyl <- as(phyl.phylo, "phylo4")
  }
  
  phyl <- reorder(phyl)
  phyl.phylo <- as(phyl, "phylo")  
  
  res <- list(phyl = phyl, phyl.phylo = phyl.phylo)
  
  }
