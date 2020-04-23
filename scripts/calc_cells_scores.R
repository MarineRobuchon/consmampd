###############################################################################################
# script to calculate spatial scores for each cell of the reference gridded map
###############################################################################################

### load packages and functions ###----
library(raster)
library(ape)
library(phytools)
library(phylobase)
library(picante)
source("./scripts/calc_PE.R")
source("./scripts/new_functions.R")
library(foreach)
library(doParallel)

### load data ###----
## scores and traits for mammals
mammals <- read.csv2(paste0(getwd(), "/outputs/mammals25092019.csv"))[, -1]

## the geoTIFF spatial files for each species in PHYLACINE
list_raster_mammals <- list.files(path = paste0(getwd(), "/requireddata/Current/Current/"), pattern = ".tif", full.names = TRUE)
stack_raster_mammals <- stack(list_raster_mammals) # takes several minutes

### calculate scores by cell ### ----
## 1. species richness
richness <- stackApply(x = subset(stack_raster_mammals, mammals$species_phylacine),
                        indices = nlayers(subset(stack_raster_mammals, mammals$species_phylacine)),
                        fun = sum) # takes about 1h
writeRaster(richness, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_richness.grd")) # save the raster
# richness <- raster(paste0(getwd(), "/outputs/Map_mammals_richness.grd"))# alternatively, directly open the calculated raster
plot(richness) # quick plot, see script figures_tables for the nice corresponding figure

## 2. threatened species richness
threatened <- mammals[which(mammals$IUCN.Status.1.2 %in% c("CR", "EN", "VU")),]
richness_threatened <- stackApply(x = subset(stack_raster_mammals, levels(factor(threatened$species_phylacine))),
                                indices = nlayers(subset(stack_raster_mammals, levels(threatened$species_phylacine))),
                                fun = sum) # takes about 15 min
writeRaster(richness_threatened, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_threatened-richness.grd"))
# richness_threatened <- raster(paste0(getwd(), "/outputs/Map_mammals_threatened-richness.grd")) # alternatively, directly open the calculated raster
plot(richness_threatened) # quick plot, see script figures_tables for the nice corresponding figure

## 3. rare species richness 
# with rare defined as binary 
# depending on the EOO of the species compared to the mediane EOO
# as in Stein et al. 2018 referring to Davy & Davidson 2017
medianEOO <- median(mammals$range_size)  # [1] 437681.3
mammals$rare <- NA
mammals$rare[mammals$range_size < medianEOO] <- "yes"
mammals$rare[mammals$range_size >= medianEOO] <- "no"
rare <- mammals[which(mammals$rare == "yes"),]
richness_rare <- stackApply(x = subset(stack_raster_mammals, levels(factor(rare$species_phylacine))),
                                 indices = nlayers(subset(stack_raster_mammals, levels(rare$species_phylacine))),
                                 fun = sum) # takes about 30 min
writeRaster(richness_rare, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_rare-richness.grd"))
# richness_rare <- raster(paste0(getwd(), "/outputs/Map_mammals_rare-richness.grd")) # alternatively, directly open the calculated raster
plot(richness_rare) # quick plot, see script figures_tables for the nice corresponding figure

## 4. species-weighted rarity
# first create a dataframe of species presence from the stack
gridded.presab <- as.data.frame(subset(stack_raster_mammals, mammals$species_phylacine)) # takes several hours
saveRDS(gridded.presab, paste0(getwd(), "/outputs/gridded.presab"))
# gridded.presab <- readRDS(paste0(getwd(), "/outputs/gridded.presab")) # alternatively, you can directly open the already generated gridded.presab
# second caculate weights of rarity for each species (1/number of cells where the species is present)
occur <- colSums(gridded.presab)
weights <- 1 / occur
# third sum weights of rarity over the species to get species-weighted rarity - takes about 1h30
weighted_rarity <- calc(subset(stack_raster_mammals, mammals$species_phylacine),
     function(x, wght = weights, ...)
     {
       sum(x * wght, na.rm = TRUE)
     })
writeRaster(weighted_rarity, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_weighted-rarity.grd"))
# weighted_rarity <- raster(paste0(getwd(), "/outputs/Map_mammals_weighted-rarity.grd")) # alternatively, directly open the calculated raster
plot(weighted_rarity) # quick plot, see script figures_tables for the nice corresponding figure

## 5. phylogenetic diversity
mam_trees <- read.nexus(paste0(getwd(), "/requireddata/Complete_phylogeny.nex"))
mam_sub100trees <- mam_trees[seq(1, length(mam_trees), 10)] # extract 100 phylogenetic trees

list_PD <- list()

for (i in 1:length(mam_sub100trees)) # takes several days to calculate PD for each grid cell and each of the 100 trees
  {
  gpd <-  picante::pd(gridded.presab, mam_sub100trees[[i]])
  gpd$PD
  phylogenetic_diversity <- subset(stack_raster_mammals, 1)
  phylogenetic_diversity[] <- NA
  phylogenetic_diversity[] <- gpd$PD
  list_PD[[i]] <- phylogenetic_diversity
  print(i)
}

stack_PD <- stack(list_PD)
median_PD <- calc(stack_PD, fun = median) # calculate the median over the 100 trees
writeRaster(median_PD, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_median-PD.grd"))
# median_PD <- raster(paste0(getwd(), "/outputs/Map_mammals_median-PD.grd")) # alternatively, directly open the calculated raster
plot(median_PD) # quick plot, see script figures_tables for the nice corresponding figure

## 6. threatened phylogenetic diversity
gridded.presab.threatened <- as.data.frame(subset(stack_raster_mammals, as.character(threatened$species_phylacine))) # takes several minutes

list_threatened_PD <- list()
 
for (i in 1:length(mam_sub100trees)) # takes several hours to calculate threatened PD for each grid cell and each of the 100 trees
  {
  gpd <-  picante::pd(gridded.presab.threatened, mam_sub100trees[[i]])
  gpd$PD
  threatened_phylogenetic_diversity <- subset(stack_raster_mammals, 1)
  threatened_phylogenetic_diversity[] <- NA
  threatened_phylogenetic_diversity[] <- gpd$PD
  list_threatened_PD[[i]] <- threatened_phylogenetic_diversity
  print(i)
}

stack_threatened_PD <- stack(list_threatened_PD)
median_threatened_PD <- calc(stack_threatened_PD, fun = median) # calculate the median over the 100 trees
writeRaster(median_threatened_PD, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_median-threatened-PD.grd"))
# median_threatened_PD <- raster(paste0(getwd(), "/outputs/Map_mammals_median-threatened-PD.grd")) # alternatively, directly open the calculated raster
plot(median_threatened_PD) # quick plot, see script figures_tables for the nice corresponding figure

## 7. rare phylogenetic diversity
gridded.presab.rare <- as.data.frame(subset(stack_raster_mammals, as.character(rare$species_phylacine))) # takes several minutes

list_rare_PD <- list()

for (i in 1:length(mam_sub100trees)) # takes several hours to calculate rare PD for each grid cell and each of the 100 trees
  {
  gpd <-  picante::pd(gridded.presab.rare, mam_sub100trees[[i]])
  gpd$PD
  rare_phylogenetic_diversity <- subset(stack_raster_mammals, 1)
  rare_phylogenetic_diversity[] <- NA
  rare_phylogenetic_diversity[] <- gpd$PD
  list_rare_PD[[i]] <- rare_phylogenetic_diversity
  print(i)
}

stack_rare_PD <- stack(list_rare_PD)
median_rare_PD <- calc(stack_rare_PD, fun = median) # calculate the median over the 100 trees
writeRaster(median_rare_PD, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_median-rare-PD.grd"))
# median_rare_PD <- raster(paste0(getwd(), "/outputs/Map_mammals_median-rare-PD.grd")) # alternatively, directly open the calculated raster
plot(median_threatened_PD) # quick plot, see script figures_tables for the nice corresponding figure

## 8. phylogenetic-weighted rarity
matrix.presab <- as.matrix(subset(stack_raster_mammals, mammals$species_phylacine))

list_PWR <- list()
 
for (i in 1:length(mam_sub100trees)) # takes several days to calculate PWR for each grid cell and each of the 100 trees
  {
  PWR <-  calc_PE(mam_sub100trees[[i]], matrix.presab, "presence")
  phylogenetic_weighted_rarity <- subset(stack_raster_mammals, 1)
  phylogenetic_weighted_rarity[] <- NA
  phylogenetic_weighted_rarity[] <- PWR$PE
  list_PWR[[i]] <- phylogenetic_weighted_rarity
  print(i)
}

stack_PWR <- stack(list_PWR)
median_PWR <- calc(stack_PWR, fun = median) # calculate the median over the 100 trees
writeRaster(median_PWR, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_median-PWR.grd"))
# median_PWR <- raster(paste0(getwd(), "/outputs/Map_mammals_median-PWR.grd")) # alternatively, directly open the calculated raster
plot(median_PWR) # quick plot, see script figures_tables for the nice corresponding figure

## 9. number of species in the TOP 25% HEDGE
tophedge <- read.csv2(paste0(getwd(), "/outputs/mammalsTOP25percenthedge.csv"))
richness_TOPHEDGE <- stackApply(x = subset(stack_raster_mammals, levels(factor(tophedge$species_phylacine))),
                                indices = nlayers(subset(stack_raster_mammals, levels(tophedge$species_phylacine))),
                                fun = sum) # takes several minutes
writeRaster(richness_TOPHEDGE, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_richnessTOPHEDGE.grd"))
# richness_TOPHEDGE <- raster(paste0(getwd(), "/outputs/Map_mammals_richnessTOPHEDGE.grd")) # alternatively, directly open the calculated raster
plot(richness_TOPHEDGE) # quick plot, see script figures_tables for the nice corresponding figure
# check which TOP HEDGE species are in/just around the Caspian Sea
stack_TOPHEDGE <- subset(stack_raster_mammals, levels(factor(tophedge$species_phylacine)))
caspian_TOPHEDGE <- crop(stack_TOPHEDGE, c(4300000, 5400000, 4300000, 5400000)) 
occ_caspian_TOPHEDGE <- cellStats(caspian_TOPHEDGE, "sum")
names(occ_caspian_TOPHEDGE[occ_caspian_TOPHEDGE!=0])

## 10. number of species in the TOP 25% LEDGE
topledge <- read.csv2(paste0(getwd(), "/outputs/mammalsTOP25percentledge.csv"))
richness_TOPLEDGE <- stackApply(x = subset(stack_raster_mammals, levels(factor(topledge$species_phylacine))),
                                indices = nlayers(subset(stack_raster_mammals, levels(topledge$species_phylacine))),
                                fun = sum) # takes several minutes
writeRaster(richness_TOPLEDGE, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_richnessTOPLEDGE.grd"))
# richness_TOPLEDGE <- raster(paste0(getwd(), "/outputs/Map_mammals_richnessTOPLEDGE.grd")) # alternatively, directly open the calculated raster
plot(richness_TOPLEDGE) # quick plot, see script figures_tables for the nice corresponding figure

## 11. expected gain in PD if all species present in a grid cell were secured (p = 0)
ref_raster <- subset(stack_raster_mammals, 1)
ref_raster[] <- NA
mam_traits <- read.csv2(paste0(getwd(), "/outputs/Trait_data_status_imputed.csv"))
mam_status <- mam_traits$Imputed.Status
mam_pext <- status2pext(mam_status)$pext_IUCN50 
names(mam_pext) <- mam_traits$Binomial.1.2
# the following lines are not run by default because the calculations take about one week
# you can either uncomment them and reproduce the calculations 
# or go directly to line 208 to load the result

# # start to calculate with parallelisation
# 
# registerDoParallel(cores = 2)
# 
# trees <- foreach(i = 1:length(mam_sub100trees),
#                  .packages = c("ape", "phytools", "phylobase", "raster"),
#                  .export = c("mam_sub100trees", "gridded.presab", "mam_pext")) %dopar%
#                  {
#                    GexpPD_p0 <-  GLexpPD.gridded (mam_sub100trees[[i]], gridded.presab, mam_pext, 0)
#                    gain_expPD_p0 <- ref_raster
#                    gain_expPD_p0[] <- GexpPD_p0
#                    saveRDS(object = gain_expPD_p0, file = paste0(getwd(), "/outputs/GexpPD_p0_tree", i))
#                  }

# # load result files, stack them, and compute the median
# list_GexpPD_p0 <- list.files(path = paste0(getwd(), "/outputs/"), pattern = "GexpPD_p0_tree", full.names = TRUE)
# rasters_GexpPD_p0 <- lapply(list_GexpPD_p0, readRDS)
# stack_GexpPD_p0 <- stack(rasters_GexpPD_p0)
# median_GexpPD_p0 <- calc(stack_GexpPD_p0, fun = median)
# writeRaster(median_GexpPD_p0, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_median-GexpPD-p0.grd"))
median_GexpPD_p0 <- raster(paste0(getwd(), "/outputs/Map_mammals_median-GexpPD-p0.grd"))
plot(median_GexpPD_p0) # quick plot, see script figures_tables for the nice corresponding figure

# identify priority areas (= the 2.5% cells with the highest expected gain in PD if all species present in a grid cell were secured)
GexpPD_hotspots <- median_GexpPD_p0
GexpPD_hotspots[GexpPD_hotspots < quantile(median_GexpPD_p0, 0.975)] <- NA # NA values for cells which are not hotspots
GexpPD_hotspots[GexpPD_hotspots >= quantile(median_GexpPD_p0, 0.975)] <- 1 # assign value 1 to priority areas
plot(GexpPD_hotspots) # quick plot, see script figures_tables for the nice corresponding figure
writeRaster(GexpPD_hotspots, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_GexpPD_hotspots.grd"))

## 12. expected loss in PD if all species present in a grid cell became extinct (p = 1)
ref_raster <- subset(stack_raster_mammals, 1)
ref_raster[] <- NA
# the following lines are not run by default because the calculations take about one week
# you can either uncomment them and reproduce the calculations 
# or go directly to line 249 to load the result

# # start to calculate with parallelisation
# rm(trees)
# registerDoParallel(cores = 2)
# trees <- foreach(i = 1:length(mam_sub100trees),
#                  .packages = c("ape", "phytools", "phylobase", "raster"),
#                  .export = c("mam_sub100trees", "gridded.presab", "mam_pext")) %dopar%
#                  {
#                    LexpPD_p1 <-  GLexpPD.gridded (mam_sub100trees[[i]], gridded.presab, mam_pext, 1)
#                    loss_expPD_p1 <- ref_raster
#                    loss_expPD_p1[] <- LexpPD_p1
#                    saveRDS(object = loss_expPD_p1, file = paste0(getwd(), "/outputs/LexpPD_p1_tree", i))
#                  }

# # load result files, stack them, and compute the median
# list_LexpPD_p1 <- list.files(path = paste0(getwd(), "/outputs/"), pattern = "LexpPD_p1_tree", full.names = TRUE)
# rasters_LexpPD_p1 <- lapply(list_LexpPD_p1, readRDS)
# stack_LexpPD_p1 <- stack(rasters_LexpPD_p1)
# median_LexpPD_p1 <- calc(stack_LexpPD_p1, fun = median)
# writeRaster(median_LexpPD_p1, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_median-LexpPD-p1.grd"))
median_LexpPD_p1 <- raster(paste0(getwd(), "/outputs/Map_mammals_median-LexpPD-p1.grd"))
plot(abs(median_LexpPD_p1)) # quick plot, see script figures_tables for the nice corresponding figure

# identify loss-significant areas (= the 2.5% cells with the highest expected loss in PD if all species present in a grid cell became extinct)
LexpPD_hotspots <- abs(median_LexpPD_p1)
LexpPD_hotspots[LexpPD_hotspots < quantile(abs(median_LexpPD_p1), 0.975)] <- NA # NA values for cells which are not hotspots
LexpPD_hotspots[LexpPD_hotspots >= quantile(abs(median_LexpPD_p1), 0.975)] <- 1 # assign value 1 to loss-significant areas
plot(LexpPD_hotspots) # quick plot, see script figures_tables for the nice corresponding figure
writeRaster(LexpPD_hotspots, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_LexpPD_hotspots.grd"))

## 13. proportion of threatened species
prop_richness_threatened <- richness_threatened/richness
plot(prop_richness_threatened) # quick plot, see script figures_tables for the nice corresponding figure
writeRaster(prop_richness_threatened, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_prop-richness-threatened.grd"))

## 14. proportion of rare species
prop_richness_rare <- richness_rare/richness
plot(prop_richness_rare) # quick plot, see script figures_tables for the nice corresponding figure
writeRaster(prop_richness_rare, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_prop-richness-rare.grd"))

## 15. mean species-weighted rarity
mean_weighted_rarity <- weighted_rarity/richness
plot(mean_weighted_rarity) # quick plot, see script figures_tables for the nice corresponding figure
writeRaster(mean_weighted_rarity, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_mean-weighted-rarity.grd"))

## 16. proportion of threatened PD
prop_threatened_PD <- median_threatened_PD/median_PD
plot(prop_threatened_PD) # quick plot, see script figures_tables for the nice corresponding figure
writeRaster(prop_threatened_PD, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_prop-threatened-PD.grd"))

## 17. proportion of rare PD 
prop_rare_PD <- median_rare_PD/median_PD
plot(prop_rare_PD) # quick plot, see script figures_tables for the nice corresponding figure
writeRaster(prop_rare_PD, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_prop-rare-PD.grd"))

## 18. mean phylogenetic-weighted rarity
mean_PWR <- median_PWR/median_PD
plot(mean_PWR) # quick plot, see script figures_tables for the nice corresponding figure
writeRaster(mean_PWR, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_mean-PWR.grd"))

## 19. proportion of TOP 25% HEDGE species
prop_richness_TOPHEDGE <- richness_TOPHEDGE/richness
plot(prop_richness_TOPHEDGE) # quick plot, see script figures_tables for the nice corresponding figure
writeRaster(prop_richness_TOPHEDGE, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_prop-richness-TOPHEDGE.grd"))

## 20. proportion of TOP 25% LEDGE species 
prop_richness_TOPLEDGE <- richness_TOPLEDGE/richness
plot(prop_richness_TOPLEDGE) # quick plot, see script figures_tables for the nice corresponding figure
writeRaster(prop_richness_TOPLEDGE, overwrite = TRUE, filename = paste0(getwd(), "/outputs/Map_mammals_prop-richness-TOPLEDGE.grd"))

