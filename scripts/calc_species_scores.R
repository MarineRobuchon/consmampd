###############################################################################################
# script to calculate HEDGE, LEDGE and ED scores
###############################################################################################

### load packages, functions and rltoken ###----
library(ape)
library(phylobase)
library(raster)
library(rgdal)
library(redlistr)
source("./scripts/new_functions.R")
source("./scripts/checkphyloarg.R")
library(adiv)
library(rredlist)
rltoken <- "d2db90179206dcdc872af96c3ec4ebba3e80d0d2e9a37bbf8831a8641a0456b1" # this is my personal token to use the IUCN Red List API, please generate your own (here: https://apiv3.iucnredlist.org/api/v3/token) and replace it
library(dplyr)

### load phylogenetic trees ###----
mam_trees <- read.nexus(paste0(getwd(), "/requireddata/Complete_phylogeny.nex"))
mam_trees <- lapply (mam_trees, FUN = function (x) {as(x, "phylo4")}) # takes several minutes
str(mam_trees)

### load IUCN status ###----
mam_traits <- read.csv(paste0(getwd(), "/requireddata/Trait_data.csv"))
str(mam_traits)

### remove species with IUCN status "EX", "EP", "EW" in mam_traits-----
mam_traits <- mam_traits[-which(mam_traits$IUCN.Status.1.2 %in% c("EX", "EP", "EW")) ,]
mam_traits$IUCN.Status.1.2 <- factor(mam_traits$IUCN.Status.1.2)
levels(mam_traits$IUCN.Status.1.2)
nrow(mam_traits) # the 5477 extant species
nrow(mam_traits[mam_traits$IUCN.Status.1.2 == "DD" ,]) # among which 777 species are DD

### imputation of status for DD species (following method of Veron et al. 2016)----
## you can either uncomment the next lines and run them to re-calculate species ranges
## or directly read "mammals_range_table.csv" with the calculated species ranges
# # calculate range size from the geoTIFF spatial files for each species in PHYLACINE
# list_raster_mammals <- list.files(path = paste0(getwd(), "/requireddata/Current/Current/"), pattern = ".tif", full.names = TRUE) # takes several minutes
# stack_raster_mammals <- stack(list_raster_mammals) # takes several minutes
# range_table <- data.frame(species = as.character(names(stack_raster_mammals)), range_size = NA)
# 
# for (i in 1:nlayers(stack_raster_mammals))
#   {
#     ras <- stack_raster_mammals[[i]]
#     range_table$species[i] <- names(ras)
#     ras.bin <- ras == 1
#     values(ras.bin)[values(ras.bin)!= 1] <- NA
#     range_table$range_size[i] <- getArea(ras.bin) # in km2
#     print(names(ras))
# }
# 
# write.csv2(range_table, paste0(getwd(), "/outputs/mammals_range_table.csv"))
range_table <- read.csv2(paste0(getwd(), "/outputs/mammals_range_table.csv"))

# for each DD species, find the 3 species with the 6 closest ranges (3 below, 3 above)
# calculate the median of the pext of these 6 species and attribute this pext to the DD species
# find the status for which the pext is the closest and attribute it to the DD species
DDspecies <- data.frame(species = as.character(mam_traits$Binomial.1.2[mam_traits$IUCN.Status.1.2 == "DD"]), att_pext = NA, att_status = NA)
DDspecies$species <- as.character(DDspecies$species)
range_table$species <- as.character(range_table$species)
range_table <- range_table[range_table$species %in% mam_traits$Binomial.1.2 ,] # we remove species which are not in mam_traits anymore (ie species EP, EX, EW)
mooers <- read.csv2(paste0(getwd(), "/data/pext_Mooers2008.csv"))

for (i in 1:nrow(DDspecies))
  {
    species <- DDspecies$species[i]
    range_table_withoutDD <- range_table[!range_table$species %in% DDspecies$species ,]
    range_table_withsp <- rbind(range_table_withoutDD, range_table[range_table$species == species ,])
    range_table_withsp <- range_table_withsp[order(range_table_withsp$range_size),]
    j <- which(range_table_withsp$species == species)
    list.species <- range_table_withsp$species[c(j+1, j+2, j+3, j-1, j-2, j-3)] 
    status.species <- mam_traits[mam_traits$Binomial.1.2 %in% list.species , c("Binomial.1.2", "IUCN.Status.1.2")] 
    pext.species <- merge(status.species, mooers, by.x = "IUCN.Status.1.2", by.y = "IUCN.Category")
    DDspecies$att_pext[i] <- median(pext.species$IUCN50)
    DDspecies$att_status[i] <- as.character(mooers$IUCN.Category[which.min(abs(mooers$IUCN50-DDspecies$att_pext[i]))])
    print(species)
  }

# now put a new colum in mam_traits with imputed status
mam_traits$Imputed.Status <- mam_traits$IUCN.Status.1.2

for (i in 1:nrow(mam_traits))
  {
    if (mam_traits$Imputed.Status[i] == "DD")
      mam_traits$Imputed.Status[i] <- DDspecies$att_status[DDspecies$species == mam_traits$Binomial.1.2[i]]
  }
write.csv2(mam_traits, paste0(getwd(), "/outputs/Trait_data_status_imputed.csv")) # safeguard imputed status together with other traits

### calculate hedge and ledge (and hed, not used in the paper) as well as ED (and other scores of phylogenetic originality, not used in the paper) ###----
## hedge, ledge and hed scores
mam_status <- mam_traits$Imputed.Status
mam_pext <- status2pext(mam_status)$pext_IUCN50 
names(mam_pext) <- mam_traits$Binomial.1.2
mam_sub100trees <- mam_trees[seq(1, length(mam_trees), 10)] # the 100 phyl inputs for the HED2 function

# the following calculations take about 5 hours
mam_HED <- list()
for (i in 1:length(mam_sub100trees))
{
  mam_tree <- mam_sub100trees[[i]]
  mam_HED[[i]] <- HED2(phyl = mam_tree, proba = mam_pext)
}

mam.HEDallscores <- lapply(X = mam_HED, function(m) m$scores) # extraction of the list of scores
mam.HEDallscores <- do.call(cbind, mam.HEDallscores) # and transform it to a dataframe
mam.hedgemed <- apply(X = mam.HEDallscores[colnames(mam.HEDallscores)=="HEDGE"], MARGIN = 1, FUN = median)  
mam.ledgemed <- apply(X = mam.HEDallscores[colnames(mam.HEDallscores)=="LEDGE"], MARGIN = 1, FUN = median)
mam.hedmed <- apply(X = mam.HEDallscores[colnames(mam.HEDallscores)=="HED"], MARGIN = 1, FUN = median)
mam.hedgesd <- apply(X = mam.HEDallscores[colnames(mam.HEDallscores)=="HEDGE"], MARGIN = 1, FUN = sd)  
mam.ledgesd <- apply(X = mam.HEDallscores[colnames(mam.HEDallscores)=="LEDGE"], MARGIN = 1, FUN = sd) 
mam.hedsd <- apply(X = mam.HEDallscores[colnames(mam.HEDallscores)=="HED"], MARGIN = 1, FUN = sd)
mam.med <- data.frame(cbind(hedge_median = mam.hedgemed, hedge_sd = mam.hedgesd, 
                            ledge_median = mam.ledgemed, ledge_sd = mam.ledgesd,
                            hed_median = mam.hedmed, hed_sd = mam.hedsd))
mam.med <- mam.med[order(mam.med$hedge_median, decreasing = TRUE),]
write.csv2(mam.med, paste0(getwd(), "/outputs/scores_Mammals_median-over-the-100-resolved-trees.csv")) # save results

## calculation of ED (as well as other measures of phylogenetic originality at the same time)
# the calculations take several hours
mam_phylori <- list()
for (i in 1:length(mam_sub100trees))
{
  mam_tree <- mam_sub100trees[[i]] 
  mam.dTree <- distinctTree(mam_tree) # ED and ES
  mam.dUltra <- distinctUltra(mam_tree) # QEbased and 2Hbased
  mam.termBL <- edgeLength(mam_tree)[getEdge(mam_tree, nodeId(mam_tree, type = "tip"))] # terminal branch length
  mam.meanD <- apply(X = cophenetic(as(mam_tree, "phylo"))/2, MARGIN = 1, FUN = mean) # mean distance to other species in the tree
  mam_phylori[[i]] <- cbind (mam.dTree, mam.dUltra, termBL = mam.termBL, meanD = mam.meanD)
}

mam.allphylori <- do.call(cbind, mam_phylori)
mam.EDmed <- apply(X = mam.allphylori[colnames(mam.allphylori)=="ED"], MARGIN = 1, FUN = median)  
mam.ESmed <- apply(X = mam.allphylori[colnames(mam.allphylori)=="ES"], MARGIN = 1, FUN = median)  
mam.QEBASEDmed <- apply(X = mam.allphylori[colnames(mam.allphylori)=="QEbased"], MARGIN = 1, FUN = median)  
mam.2HBASEDmed <- apply(X = mam.allphylori[colnames(mam.allphylori)=="2Hbased"], MARGIN = 1, FUN = median)  
mam.termBLmed <- apply(X = mam.allphylori[colnames(mam.allphylori)=="termBL"], MARGIN = 1, FUN = median)  
mam.meanDmed <- apply(X = mam.allphylori[colnames(mam.allphylori)=="meanD"], MARGIN = 1, FUN = median) 
mam.EDsd <- apply(X = mam.allphylori[colnames(mam.allphylori)=="ED"], MARGIN = 1, FUN = sd)  
mam.ESsd <- apply(X = mam.allphylori[colnames(mam.allphylori)=="ES"], MARGIN = 1, FUN = sd)  
mam.QEBASEDsd <- apply(X = mam.allphylori[colnames(mam.allphylori)=="QEbased"], MARGIN = 1, FUN = sd)  
mam.2HBASEDsd <- apply(X = mam.allphylori[colnames(mam.allphylori)=="2Hbased"], MARGIN = 1, FUN = sd)  
mam.termBLsd <- apply(X = mam.allphylori[colnames(mam.allphylori)=="termBL"], MARGIN = 1, FUN = sd)  
mam.meanDsd <- apply(X = mam.allphylori[colnames(mam.allphylori)=="meanD"], MARGIN = 1, FUN = sd) 
mam_PHYLORImed <- data.frame(cbind(median_ED = mam.EDmed, sd_ED = mam.EDsd,
                                       median_ES = mam.ESmed, sd_ES = mam.ESsd,
                                       median_QEbased = mam.QEBASEDmed, sd_QEbased = mam.QEBASEDsd,
                                       median_2Hbased = mam.2HBASEDmed, sd_2Hbased = mam.2HBASEDsd,
                                       median_termBL = mam.termBLmed, sd_termBL = mam.termBLsd,
                                       median_meanD = mam.meanDmed, sd_meanD = mam.meanDsd ))
write.csv2(mam_PHYLORImed, paste0(getwd(), "/outputs/phylori_Mammals_median-over-the-100-resolved-trees.csv")) # save results

### create a dataframe merging all calculated species scores and species traits ----
### and distinguish between names used in Phylacine (IUCN 2016-3) and names accepted in the current IUCN version (2019-2 at the time of calculations)
# mam.med <- read.csv2(paste0(getwd(), "/outputs/scores_Mammals_median-over-the-100-resolved-trees.csv")) # run if you did not manage to generate mam.med yourself
mam.med <- mam.med[, c(1:2, 4, 6)] # select the median scores
colnames(mam.med) <- c("species", "hedge", "ledge", "hed") # rename the columns
mam.med$hedge_rank <- rank(-mam.med$hedge, ties.method = "average") # calculate ranks
mam.med$ledge_rank <- rank(-mam.med$ledge, ties.method = "average")
mam.med$hed_rank <- rank(-mam.med$hed, ties.method = "average")

# mam_PHYLORImed <- read.csv2(paste0(getwd(), "/outputs/phylori_Mammals_median-over-the-100-resolved-trees.csv")) # run if you did not manage to generate mam_PHYLORImed yourself
mam_PHYLORImed <- mam_PHYLORImed[, c(1:2, 4, 6, 8, 10, 12)] # select the median scores
colnames(mam_PHYLORImed) <- c("species", "ED", "ES", "QEbased", "2Hbased", "termBL", "meanD") # rename the columns

## merge scores and traits in a same dataframe
mammals <- merge(mam.med, mam_PHYLORImed, by = "species")
mammals <- merge(mammals, mam_traits, by.x = "species", by.y = "Binomial.1.2")
mammals <- merge(mammals, range_table[, -1], by = "species")

## distinguish between names used in Phylacine (IUCN 2016-3) and names accepted in the current IUCN version (2019-2 at the time of calculations)
rl_version(key = rltoken) # current IUCN version
colnames(mammals)[1] <- "species_phylacine"
mammals_currentIUCN <- read.csv(paste0(getwd(), "/data/assessments_mammals_25092019.csv")) # reference for the 2019-2 version
mammals_currentIUCN$scientificName # list of accepted names in the 2019-2 version

mammals$species_currentIUCN <- NA

for (i in 1:nrow(mammals)) # this takes several minutes
{
  species2search <-  as.character(gsub(pattern = "_", replacement = " ", x = mammals$species_phylacine[i]))
  if (species2search %in% mammals_currentIUCN$scientificName)
  {
    result <- NA
    mammals$species_currentIUCN[i] <- species2search
  }
  else
  {
    result <- rl_synonyms(name = species2search, key = rltoken)$result
    if (nlevels(factor(result$accepted_id)) == 1)
    {
    mammals$species_currentIUCN[i] <- result$accepted_name[1]
    }
    else
    {
    mammals$species_currentIUCN[i] <- "no/not easy equivalence with phylacine"
    }
  }
}

mammals$species_phylacine[which(mammals$species_currentIUCN == "no/not easy equivalence with phylacine")] # 9 problematic names
rl_synonyms(name = "Callicebus_cupreus", key = rltoken)
mammals$species_currentIUCN[which(mammals$species_phylacine == "Callicebus_cupreus")] <- "Plecturocebus cupreus"
rl_synonyms(name = "Callicebus_torquatus", key = rltoken)
mammals$species_currentIUCN[which(mammals$species_phylacine == "Callicebus_torquatus")] <- "Cheracebus torquatus"
rl_synonyms(name = "Cebus_capucinus", key = rltoken) # on IUCN website, 2 subspecies without synonyms - we leave it uncorrected
rl_synonyms(name = "Chelemys_delfini", key = rltoken) # on MSW3, needs taxonomic refinement:
# some authors considered it's a subspecies of C. meaglonys, a name already taken in our list of species
# so we leave it uncorrected
rl_synonyms(name = "Mormopterus_loriae", key = rltoken)
mammals$species_currentIUCN[which(mammals$species_phylacine == "Mormopterus_loriae")] <- "Ozimops loriae"
rl_synonyms(name = "Mus_orangiae", key = rltoken) # # on MSW3, needs taxonomic refinement:
# Either treated as a species possibly allied to M. setzeri or listed as a subspecies of M. minutoides
# both names already taken so we leave it uncorrected
rl_synonyms(name = "Otomys_saundersiae", key = rltoken)
mammals$species_currentIUCN[which(mammals$species_phylacine == "Otomys_saundersiae")] <- "Otomys karoensis"
rl_synonyms(name = "Rattus_arfakienis", key = rltoken)
mammals$species_currentIUCN[which(mammals$species_phylacine == "Rattus_arfakienis")] <-"Rattus arfakiensis"
rl_synonyms(name = "Sorex_arunchi", key = rltoken) # on MSW3, needs taxonomic refinement:
# Possibly related or conspecific with Sorex antinorii
# name already taken so we leave it uncorrected
write.csv2(mammals, paste0(getwd(), "/outputs/mammals25092019.csv")) # save the final file

## all the previous steps of data merging and names equivalence (lines 160 to 227) can be replaced by the following line
mammals <- read.csv2(paste0(getwd(), "/outputs/mammals25092019.csv"))[, -1]

### identify species of interest ### ----
## 1. identification of the TOP 25%  HEDGE species
tophedge <- head(mammals[order(mammals$hedge_rank) ,], n = round(nrow(mammals)*25/100))
head(tophedge)
write.csv2(tophedge, paste0(getwd(), "/outputs/mammalsTOP25percenthedge.csv"))

## 2. identification of the TOP 25% LEDGE species
topledge <- head(mammals[order(mammals$ledge_rank) ,], n = round(nrow(mammals)*25/100))
head(topledge)
write.csv2(topledge, paste0(getwd(), "/outputs/mammalsTOP25percentledge.csv"))

### conservation measures and introduction status of species of interest ### ----
## 1. for the TOP 25%  HEDGE species
# first list conservation measures
result_list <- list()

for (i in 1:nrow(tophedge)) # takes about 100 min
{
  species <-  as.character(tophedge$species_currentIUCN[i])
  result <- rl_measures(name = species, key = rltoken)$result
  result$species_phylacine <- tophedge$species_phylacine[i]
  result$species_currentIUCN <- tophedge$species_currentIUCN[i]
  result_list[[i]] <- result
  print(i)
}

result_tophedge <- bind_rows(result_list)
alldata_tophedge <- merge(tophedge, result_tophedge, by = "species_phylacine")
alldata_tophedge$code[is.na(alldata_tophedge$code)] <- "none"
alldata_tophedge$title[is.na(alldata_tophedge$title)] <- "none"
levels(factor(alldata_tophedge$code))
levels(factor(alldata_tophedge$title))
alldata_tophedge$conservation_classification <-  NA
alldata_tophedge$conservation_classification[substr(alldata_tophedge$code, 1, 2)=="1."] <- "1 Land/water protection"
alldata_tophedge$conservation_classification[substr(alldata_tophedge$code, 1, 2)=="2."] <- "2 Land/water management"
alldata_tophedge$conservation_classification[substr(alldata_tophedge$code, 1, 2)=="3."] <- "3 Species management"
alldata_tophedge$conservation_classification[substr(alldata_tophedge$code, 1, 2)=="4."] <- "4 Education & awareness"
alldata_tophedge$conservation_classification[substr(alldata_tophedge$code, 1, 2)=="5."] <- "5 Law & policy"
alldata_tophedge$conservation_classification[substr(alldata_tophedge$code, 1, 2)=="6."] <- "6 Livelihood, economic & other incentives"
alldata_tophedge$conservation_classification[alldata_tophedge$code == "none"] <- "None"
traits_imputed <- read.csv2(paste0(getwd(), "/outputs/Trait_data_status_imputed.csv"))
alldata_tophedge <- merge(alldata_tophedge, traits_imputed[, c("Binomial.1.2", "Imputed.Status")], by.x = "species_phylacine", by.y = "Binomial.1.2")

# second assign introduction status
# load the distribution data from IUCN for mammals (2019-2 version) 
shp_iucn_mammals <- shapefile(paste0(getwd(), "/requireddata/MAMMALS.shp")) # takes several minutes
shp_iucn_tophedge <- shp_iucn_mammals[shp_iucn_mammals@data$binomial %in% alldata_tophedge$species_currentIUCN.x ,]
levels(factor(shp_iucn_tophedge@data$origin))
# so to have a list of species introduced at least once, we select the names of polygons having "3" as origin
shp_iucn_tophedge_introduced <- shp_iucn_tophedge[shp_iucn_tophedge@data$origin == "3" ,]
tophedge_introduced <- levels(factor(shp_iucn_tophedge_introduced@data$binomial))
alldata_tophedge$introduced <- "no"
alldata_tophedge$introduced[alldata_tophedge$species_currentIUCN.x %in% tophedge_introduced] <- "yes"
write.csv2(alldata_tophedge, paste0(getwd(), "/outputs/mammals_tophedge_conservation_introduction.csv"))

## 2. for the TOP 25%  LEDGE species
# first list conservation measures
result_list <- list()

for (i in 1:nrow(topledge)) # takes about 100 min
{
  species <-  as.character(topledge$species_currentIUCN[i])
  if (species == "no/not easy equivalence with phylacine")
    result <- NA
  else
  result <- rl_measures(name = species, key = rltoken)$result
  result$species_phylacine <- topledge$species_phylacine[i]
  result$species_currentIUCN <- topledge$species_currentIUCN[i]
  result_list[[i]] <- result
  print(i)
}
 
result_topledge <- bind_rows(result_list)
alldata_topledge <- merge(topledge, result_topledge, by = "species_phylacine")
alldata_topledge$code[is.na(alldata_topledge$code)] <- "none"
alldata_topledge$title[is.na(alldata_topledge$title)] <- "none"
levels(factor(alldata_topledge$code))
levels(factor(alldata_topledge$title))
alldata_topledge$conservation_classification <-  NA
alldata_topledge$conservation_classification[substr(alldata_topledge$code, 1, 2)=="1."] <- "1 Land/water protection"
alldata_topledge$conservation_classification[substr(alldata_topledge$code, 1, 2)=="2."] <- "2 Land/water management"
alldata_topledge$conservation_classification[substr(alldata_topledge$code, 1, 2)=="3."] <- "3 Species management"
alldata_topledge$conservation_classification[substr(alldata_topledge$code, 1, 2)=="4."] <- "4 Education & awareness"
alldata_topledge$conservation_classification[substr(alldata_topledge$code, 1, 2)=="5."] <- "5 Law & policy"
alldata_topledge$conservation_classification[substr(alldata_topledge$code, 1, 2)=="6."] <- "6 Livelihood, economic & other incentives"
alldata_topledge$conservation_classification[alldata_topledge$code == "none"] <- "None"
alldata_topledge <- merge(alldata_topledge, traits_imputed[, c("Binomial.1.2", "Imputed.Status")], by.x = "species_phylacine", by.y = "Binomial.1.2")

# second assign introduction status
shp_iucn_topledge_introduced <- shp_iucn_topledge[shp_iucn_topledge@data$origin == "3" ,]
topledge_introduced <- levels(factor(shp_iucn_topledge_introduced@data$binomial))
alldata_topledge$introduced <- "no"
alldata_topledge$introduced[alldata_topledge$species_currentIUCN.x %in% topledge_introduced] <- "yes"
levels(factor(alldata_topledge$introduced))
write.csv2(alldata_topledge, paste0(getwd(), "/outputs/mammals_topledge_conservation_introduction.csv"))
