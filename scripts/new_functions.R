###############################################################################################
# script containing all new functions created for the analyses
###############################################################################################

### function 'status2pext' to transform IUCN status in probabilities of extinction
### according to the different modalities presented in Mooers et al. 2008
# status = named vector of IUCN status for each species (names are species names)
# the result is a dataframe with column 1 corresponding to the status a
# and the following columns to pext using different modalities
status2pext <- function (status) 
{
  mooers <- read.csv2(paste0(getwd(), "/data/pext_Mooers2008.csv"), stringsAsFactors = FALSE)
  mooers$GE <- c(0, 1, 2, 3, 4)
  statusprob <- data.frame(IUCN_status = as.character(status), pext_Isaac100 = NA, pext_IUCN100 = NA, pext_IUCN50 = NA, GE = NA)
  rownames(statusprob) <- names(status)
  for (i in 1:nrow(statusprob))
  {
    if (statusprob$IUCN_status[i] %in% c("DD", "EX", "EP", "EW", ""))
    {
      statusprob$pext_Isaac100[i] <- NA
      statusprob$pext_IUCN100[i] <- NA
      statusprob$pext_IUCN50[i] <- NA
      statusprob$pext_GE[i] <- NA
    }
    else
    {
      statusprob$pext_Isaac100[i] <- mooers$Isaac100[statusprob$IUCN_status[i]==mooers$IUCN.Category]
      statusprob$pext_IUCN100[i] <- mooers$IUCN100[statusprob$IUCN_status[i]==mooers$IUCN.Category]
      statusprob$pext_IUCN50[i] <- mooers$IUCN50[statusprob$IUCN_status[i]==mooers$IUCN.Category]
      statusprob$GE[i] <- mooers$GE[statusprob$IUCN_status[i]==mooers$IUCN.Category]
    }
  }
  return(statusprob)
}


### HED2: function to calculate HEDGE, LEDGE and HED scores 
### which also returns the pruned tree used for calculating scores and the scores by node

HED2 <- function(phyl, proba, subtree=FALSE, tol=1e-8){
  
  # phyl is an object of class hclust, phylo or phylo4
  # proba is a vector with names corresponding to the tips of tree
  # the function itself put the species in the same order in the phylogenetic tree and in proba
  # If proba does not contain all tips of phyl, then phyl is pruned and a new root is defined at the most recent common ancestor of the retained species as the root of the pruned tree. If subtree is FALSE, the branch between the new root and the root of phyl is kept; if subtree is TRUE, it is removed. 
  if(is.null(names(proba)) | any(is.na(names(proba)))) stop("Vector proba must have names corresponding to the tips of the phylogenetic tree")
  nsp <- length(proba)
  if(nsp==1) stop("At least two tips should be selected in proba")
  arg.phyl <- .checkphyloarg(phyl)
  tre <- arg.phyl$phyl.phylo
  if(!is.rooted(tre)) stop("phyl must be a rooted phylogenetic tree")
  if(any(proba < -tol) | any(is.na(proba))) stop("Values in proba must be known and positive")
  if(!all(names(proba)%in%tre$tip.label)) stop("missing tips in the phylogenetic tree in comparison with names of vector proba")
  if(subtree){
    tre <- drop.tip(tre, tip=tre$tip.label[!tre$tip.label%in%names(proba)])
  } else
  {
    tre <- drop.tip(tre, tip=tre$tip.label[!tre$tip.label%in%names(proba)], root.edge = tre$Nnode)
  } 
  tre4 <- as(tre, "phylo4")
  if(!hasNodeLabels(tre4)){
    nodeLabels(tre4) <- names(nodeLabels(tre4))
  } else{
    e <- nodeLabels(tre4)
    e[is.na(e)] <- names(e[is.na(e)])
    nodeLabels(tre4) <- e
  }
  proba <- proba[tipLabels(tre4)]
  if(nNodes(tre4)==1){
    des <- descendants(tre4, nodeLabels(tre4), type="tips")
    prob <- prod(proba[names(des)])
  } else{
    des <- descendants(tre4, nodeLabels(tre4), type="tips")
    prob <- unlist(lapply(des, function(x) prod(proba[x])))
  }
  valnodes <- prob*edgeLength(tre4, nodeLabels(tre4))
  valtips <- proba*edgeLength(tre4, names(proba))
  names(valtips) <- names(proba)
  names(valnodes) <- nodeLabels(tre4)
  vals <- c(valnodes, valtips)
  # score for the tips (= the species)
  hedge <- unlist(lapply(ancestors(tre4, tipLabels(tre4), "ALL"), function(x) sum(vals[names(x)], na.rm=TRUE)))
  ledge <- hedge/proba[names(hedge)]*(1-proba[names(hedge)])
  hed <- hedge + ledge
  res <- cbind.data.frame(HEDGE=hedge, LEDGE=ledge, HED=hed)
  rownames(res) <- names(hedge)
  # scores for all the nodes (including tips)
  Nhedge <- unlist(lapply(ancestors(tre4, getNode(tre4), "ALL"), function(x) sum(vals[names(x)], na.rm=TRUE)))
  Nproba <- c(proba, prob)
  Nledge <- Nhedge/Nproba*(1-Nproba)
  Nhed <- Nhedge + Nledge
  Nres <- cbind.data.frame(HEDGE=Nhedge, LEDGE=Nledge, HED=Nhed)
  rownames(Nres) <- names(Nhedge)
  # list of outputs
  reslist <- list(tre4, res, Nres)
  names(reslist) <- c("pruned.tree", "scores", "node.scores")
  return(reslist)
}



### function 'GLexpPD.gridded' to calculate gain and losses in expected PD
# if all species in a cell become extinct or are saved from extinction
# tree = global phylogenetic tree
# gridded.presab = data frame with rows as grid cell numbers matching the rasters and columns as 'species' (layers) with 0/1 values for presence of that species in that cell
# preal = real extinction probabilities for all species in the global phylogenetic tree
# pcond = conditional extinction probabilities for each species in the regional species pool (= 0 to get gain in expPD, 1 to get loss in expPD)
# depends on another function, the function cellval

GLexpPD.gridded <- function (tree, gridded.presab, preal, pcond) {
  tre4 <- as(tree, "phylo4")
  if(!all(row.names(preal)%in%tipLabels(tre4))) 
  {
    stop("missing tips in the phylogenetic tree")
  }
  # check if there are species in preal which are not in tre4
  tre4 <- prune(tre4, tips.exclude=tipLabels(tre4)[!tipLabels(tre4)%in%names(preal)])
  # remove species in tre4 which are not in preal
  if(!hasNodeLabels(tre4))
  {
    nodeLabels(tre4) <- names(nodeLabels(tre4))
  } else {
    e <- nodeLabels(tre4)
    e[is.na(e)] <- names(e[is.na(e)])
    nodeLabels(tre4) <- e
  }
  preal <- preal[match(tipLabels(tre4), names(preal))]
  # this makes sure that species are in the same order in preal and in tre4
  des <- descendants(tre4, nodeLabels(tre4), type="tips")
  prob <- unlist(lapply(des, function(x) prod(preal[][x]))) 
  # product of extinction proba for internal branches
  
  # Long calculations:
  # internal branch length
  internal.length <- edgeLength(tre4, nodeLabels(tre4)) 
  # tip length
  tip.length <- edgeLength(tre4, tipLabels(tre4)) 
  
  # product of extinction proba for internal branches
  valnodes <- prob * internal.length
  # extinction proba by egdge length for internal branches
  valtips <- preal[] * tip.length
  # extinction probas by edge length for tips
  names(valtips) <- tipLabels(tre4)
  names(valnodes) <- nodeLabels(tre4)
  vals <- c(valnodes, valtips)
  expPDloss <- sum(vals, na.rm = TRUE)
  result <- apply(gridded.presab, 1, cellval,
                  total.tree = tre4, 
                  desc = des, 
                  probas.ext = preal, 
                  cond.probas = pcond,
                  internal.branch.length = internal.length,
                  tip.branch.length = tip.length,
                  initial.expPDloss = expPDloss)
  return(result)
}


cellval <- function(cell.presab,
                    total.tree,
                    desc,
                    probas.ext,
                    cond.probas,
                    internal.branch.length,
                    tip.branch.length,
                    initial.expPDloss) 
{
  presi <- cell.presab
  pprim <- probas.ext 
  pprim[which(names(pprim) %in% names(presi)[which(presi==1)])] <- cond.probas
  # vector of extinction proba under new conditions
  probprim <- unlist(lapply(desc, function(x) prod(pprim[][x]))) 
  # product of extinction proba for internal branches under new conditions
  valnodesprim <- probprim * internal.branch.length
  # extinction proba by egdge length for internal branches under new conditions
  valtipsprim <- pprim[] * tip.branch.length
  # extinction probas by edge length for tips under new conditions
  # names(valtipsprim) <- tipLabels(tre4)
  # names(valnodesprim) <- nodeLabels(tre4)
  valsprim <- c(valnodesprim, valtipsprim)
  expPDlossprim <- sum(valsprim, na.rm = TRUE)
  score <- initial.expPDloss - expPDlossprim
  return(score)
}

### function 'corPlot' to plot pairwise correlations between variables
### modified from the function of the same name in the R package Rarity
corPlot <- function(df, method = "spearman", digits = 2, na.action = "keep", ties.method = "average",
                    title = "", xlab = "variable.name", ylab = "variable.name", ...)
{
  rankDf <- apply(df, 2, .rankings, ties.method = ties.method, na.action = na.action)
  
  if(length(xlab) > 1) if(length(xlab) != (ncol(df))) stop("Provide as many labels for x axes as there are variables, in the same order as the column order.")
  if(length(ylab) > 1) if(length(ylab) != (ncol(df))) stop("Provide as many labels for y axes as there are variables, in the same order as the column order.")
  
  crossCongR <- NULL
  crossCongP <- NULL
  crossCongICINF <- NULL
  crossCongICSUP <- NULL
  
  for (i in seq_len(ncol(rankDf)))
  {
    r <- NULL
    p <- NULL
    icinf <- NULL
    icsup <- NULL
    for (j in 1:ncol(df))
    {
      r <- c(r, cor.test(df[, i], df[, j], method = method)$estimate)
      p <- c(p, cor.test(df[, i], df[, j], method = method)$p.value)
      if (method == "spearman")
      {
      icinf <- c(icinf, spearman.ci(df[, i], df[, j], nrep = 1000, conf.level = 0.95)$conf.int[1])
      icsup <- c(icsup, spearman.ci(df[, i], df[, j], nrep = 1000, conf.level = 0.95)$conf.int[2])
      }
      else
      {
      icinf <- c(icinf, cor.test(df[, i], df[, j], method = method, conf.level = 0.95)$conf.int[1])
      icsup <- c(icsup, cor.test(df[, i], df[, j], method = method, conf.level = 0.95)$conf.int[2])
      }
    }
    crossCongR <- rbind(crossCongR, r)
    crossCongP <- rbind(crossCongP, p)
    crossCongICINF <- rbind(crossCongICINF, icinf)
    crossCongICSUP <- rbind(crossCongICSUP, icsup)
  }
  colnames(crossCongP) <- colnames(rankDf)[1:ncol(rankDf)]
  rownames(crossCongP) <- colnames(rankDf)
  colnames(crossCongR) <- colnames(rankDf)[1:ncol(rankDf)]
  rownames(crossCongR) <- colnames(rankDf)
  colnames(crossCongICINF) <- colnames(rankDf)[1:ncol(rankDf)]
  rownames(crossCongICINF) <- colnames(rankDf)
  colnames(crossCongICSUP) <- colnames(rankDf)[1:ncol(rankDf)]
  rownames(crossCongICSUP) <- colnames(rankDf)
  
  crossCongR <- round(crossCongR, 2)
  crossCongICINF <- round(crossCongICINF, 2)
  crossCongICSUP <- round(crossCongICSUP, 2)
  
  par(mfrow = c(ncol(df), ncol(df)), oma = c(4.1, 4.1, 0.6, 0.6), mar = c(.5, .5, .5, .5))
  
  for (i in 1:ncol(df))
  {
    for (j in 1:ncol(df))
    {
      if (j == i)
      {
        plot(1, 1, pch ="", bty = "o", xaxt = "n", yaxt = "n")
        text(x = 1, y = 1, paste(ifelse(is.null(colnames(df)),
                                        paste("Var", i, sep = "."),
                                        colnames(df)[i]),
                                 "\nn = ", length(which(!is.na(df[, i])))), 
             font = 2)
      } else
        if (j > i)
        {
          plot(1, 1, pch ="", xaxt = "n", yaxt = "n", bty = "n")
          text(x = 1, y = 1, paste(formatC(crossCongR[i, j], format = "f", digits = digits), 
                                   "\ninf = ", formatC(crossCongICINF[i, j], format = "f", digits = digits),
                                   "\nsup = ", formatC(crossCongICSUP[i, j], format = "f", digits = digits),
                                   if(crossCongP[i, j] < 0.0001) {"\np < 0.0001"}
                                   else {paste0("\np = ",formatC(crossCongP[i, j], format = "f"))}), font = 2)
        } else
        {
          if (length(ylab) > 1)
          {
            if (j == 1) mtext(paste(ylab[i]), outer = T, side = 2, line = 2, at = (ncol(df)-i + .5)/ncol(df), cex = .8)
          } else if (ylab == "variable.name")
          {
            if (j == 1) mtext(paste(ifelse(is.null(colnames(rankDf)),
                                           ifelse(method == "pearson",
                                                  paste("Var. ", i, sep = ""),
                                                  paste("Var. ", i, " rank", sep = "")),
                                           ifelse(method == "pearson",
                                                  colnames(rankDf)[i],
                                                  paste(colnames(rankDf)[i], " rank", sep = "")))), outer = T, side = 2, line = 2, at = (ncol(df)-i + .5)/ncol(df), cex = .8)
          } else
          {
            if (j == 1) mtext(paste(ylab), outer = T, side = 2, line = 2, at = (ncol(df)-i + .5)/ncol(df), cex = .8)
          }
          if (length(xlab) > 1)
          {
            if (i == ncol(df)) mtext(paste(xlab[j]), outer = T, side = 1, line = 2, at = (j - .5)/ncol(df), cex = .8)
          } else if(xlab == "variable.name")
          {
            if (i == ncol(df)) mtext(paste(ifelse(is.null(colnames(rankDf)),
                                                  ifelse(method == "pearson",
                                                         paste("Var. ", j, sep = ""),
                                                         paste("Var. ", j, " rank", sep = "")),
                                                  ifelse(method == "pearson",
                                                         colnames(rankDf)[j],
                                                         paste(colnames(rankDf)[j], " rank", sep = "")))), outer = T, side = 1, line = 2, at = (j - .5)/ncol(df), cex = .8)
          } else
          {
            if (i == ncol(df)) mtext(paste(xlab), outer = T, side = 1, line = 2, at = (j - .5)/ncol(df), cex = .8)
          }
          if(method == "spearman" | method == "kendall")
          {
            plot(rankDf[, i] ~ rankDf[, j], xaxt = "n", yaxt = "n", xlim = c(0, nrow(df)), ylim = c(0, nrow(df)), xlab = "", ylab = "", ...)
            if (i == ncol(df)) axis(1, at = c(0, round(nrow(df)/2, 0), nrow(df)), las = 1)
            if (j == 1) axis(2, at = c(0, round(nrow(df)/2, 0), nrow(df)), las = 1)
          } else if (method == "pearson")
          {
            dx <- diff(c(min(df[, j], na.rm = T), max(df[, j], na.rm = T)))
            dy <- diff(c(min(df[, i], na.rm = T), max(df[, i], na.rm = T)))
            plot(df[, i] ~ df[, j], xaxt = "n", yaxt = "n", xlim = c(min(df[, j], na.rm = T) - 0.05 * dx, 
                                                                     max(df[, j], na.rm = T) + 0.05 * dx), 
                 ylim = c(min(df[, i], na.rm = T) - 0.05 * dy, 
                          max(df[, i], na.rm = T) + 0.05 * dy), 
                 xlab = "", ylab = "", ...)
            if (i == ncol(df)) axis(1, las = 1)
            if (j == 1) axis(2, las = 1)
          }
        }
    }
  }
  mtext(paste(title), outer = T, side = 3, line = -1, at = .55, cex = 1.2)
}

.rankings <- function(x, ties.method = "average", na.action = "keep") length(x) - length(which(is.na(x))) + 1 - rank(x, ties.method = ties.method, na.last = na.action)







