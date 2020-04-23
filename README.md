Re-identifying species and areas of interest to conserve global mammalian phylogenetic diversity
================
Github repository & scripts created by Marine Robuchon, with contributions of Boris Leroy & Sandrine Pavoine

-   [1. Description & organisation of the repository](#description-organisation-of-the-repository)
    -   [1.1. data](#data)
    -   [1.2. images](#images)
    -   [1.3. requireddata (only on your computer)](#requireddata-only-on-your-computer)
    -   [1.4. scripts](#scripts)
    -   [1.5. outputs](#outputs)
-   [2. Requirements](#requirements)
    -   [2.1. clone repository and create additional folders](#clone-repository-and-create-additional-folders)
    -   [2.2. download data](#download-data)
    -   [2.3. install packages](#install-packages)
-   [3. Analyses](#analyses)
    -   [3.1. step 1: calculate species scores and identify species of interest](#step-1-calculate-species-scores-and-identify-species-of-interest)
    -   [3.2. step 2: calculate cells scores and identify areas of interest](#step-2-calculate-cells-scores-and-identify-areas-of-interest)
    -   [3.3. step 3: generate figures and tables](#step-3-generate-figures-and-tables)

1. Description & organisation of the repository
-----------------------------------------------

The purpose of this github repository is to document the data used, the analyses carried out and the figures and tables produced with RStudio for the paper "Re-identifying species and areas of interest to conserve global mammalian phylogenetic diversity". The repository contains 5 folders (data, images, requireddata, scripts, outputs), which content is described below.

### 1.1. data

This folder includes the following files:

-   pext\_Mooers2008.csv (file containing the equivalence between IUCN status and probabilities of extinction according to the model of Mooers et al. 2008)
-   assessments\_mammals\_25092019.csv (file containing global assessments for mammal species in IUCN version 2019-2)
-   ne\_50m\_coastline.shp, .shx, .prj, .dbf, .cpg (spatial files for the coastline from Natural Earth <https://www.naturalearthdata.com/downloads/50m-physical-vectors/>)
-   percentage\_pa\_grid\_raster.tif (raster containing the percentage of area covered by the re-processed and transformed WDPA dataset for each grid cell of the reference gridded map).

### 1.2. images

This folder includes the following subfolders:

-   raw (folder containing the original images used in the figures, together with a table listing the source and credit for each image)
-   top10hedge (folder containing transformed images of the TOP 10 HEDGE species)
-   top10ledge (folder containing transformed images of the TOP 10 LEDGE species)
-   introduced (folder containing transformed images of some species of interest that are also introduced).

### 1.3. requireddata (only on your computer)

This folder contains data too heavy to be stored in the online github repository, so you will have to create it on your computer and download data yourself (see section "Requirements" below). Once this is done, this folder should include the following files and folders:

-   Complete\_phylogeny.nex (file containing phylogenetic trees)
-   Trait\_data.csv (file containing trait data)
-   Current (folder containing species range maps)
-   MAMMALS.shp, .shx, .sbx, .sbn, .prj, .dbf, .cpg (shapefiles of mammal species distributions from IUCN version 2019-2).

### 1.4. scripts

This folder includes the following scripts:

-   calc\_species\_scores.R (script to calculate HEDGE, LEDGE and ED scores)
-   new\_functions.R (script containing all new functions created for the analyses)
-   checkphyloarg.R (function to check the arguments of a phylo object)
-   calc\_cells\_scores.R (script to calculate all spatial scores for each cell of the reference gridded map)
-   calc\_PE.R (function to calculate phylogenetic-weighted rarity, slightly modified from <https://github.com/DanRosauer/phylospatial/tree/master/PhyloEndemism_in_R>)
-   figures\_tables.R (script to generate all figures and tables of the paper).

### 1.5. outputs

This folder contains the outputs of the analyses and you should be able to generate them yourself. They are however included here in case you do not manage to generate them yourself. They include the following files:

-   mammals\_range\_table.csv (file containing the range size calculated for mammal species)
-   Trait\_data\_status\_imputed.csv (file containing trait data plus imputed status)
-   scores\_Mammals\_median-over-the-100-resolved-trees.csv (file containing median HEDGE, LEDGE and HED scores over the 100 trees)
-   phylori\_Mammals\_median-over-the-100-resolved-trees.csv (file containing median scores of phylogenetic originality over the 100 trees)
-   mammals25092019.csv (file containing all scores and traits together with a comparison of species names between Phylacine and IUCN 2019-2)
-   mammalsTOP25percenthedge.csv (file containing all scores and traits for the TOP 25% HEDGE)
-   mammalsTOP25percentledge.csv (file containing all scores and traits for the TOP 25% LEDGE)
-   mammals\_tophedge\_conservation\_introduction.csv (file listing conservation measures and introduction status - in addition to all scores and traits - for the TOP 25% HEDGE)
-   mammals\_topledge\_conservation\_introduction.csv (file listing conservation measures and introduction status - in addition to all scores and traits - for the TOP 25% LEDGE)
-   Map\_mammals\_richness.grd, .gri (gridded map of species richness)
-   Map\_mammals\_threatened-richness.grd, .gri (gridded map of threatened species richness)
-   Map\_mammals\_rare-richness.grd, .gri (gridded map of rare species richness)
-   gridded.presab (dataframe containing the presence-absence of the species for each cell of the reference gridded map)
-   Map\_mammals\_weighted-rarity.grd, .gri (gridded map of species-weighted rarity)
-   Map\_mammals\_median-PD.grd, .gri (gridded map of median phylogenetic diversity over 100 trees)
-   Map\_mammals\_median-threatened-PD.grd, .gri (gridded map of median threatened phylogenetic diversity over 100 trees)
-   Map\_mammals\_median-rare-PD.grd, .gri (gridded map of median rare phylogenetic diversity over 100 trees)
-   Map\_mammals\_median-PWR.grd, .gri (gridded map of median phylogenetic-weighted rarity over 100 trees)
-   Map\_mammals\_richnessTOPHEDGE.grd, .gri (gridded map of the number of species in the TOP 25% HEDGE)
-   Map\_mammals\_richnessTOPLEDGE.grd, .gri (gridded map of the number of species in the TOP 25% LEDGE)
-   100 files named with the pattern GexpPD\_p0\_treei, with i going from 1 to 100 (expected gain in PD if all species present in a grid cell were secured for each of the 100 trees)
-   Map\_mammals\_median-GexpPD-p0.grd, .gri (gridded map of median expected gain in PD if all species present in a grid cell were secured over 100 trees)
-   Map\_mammals\_GexpPD\_hotspots.grd, .gri (raster of priority areas = the 2.5% cells with the highest expected gain in PD if all species present in a grid cell were secured)
-   100 files named with the pattern LexpPD\_p1\_treei, with i going from 1 to 100 (expected loss in PD if all species present in a grid cell became extinct for each of the 100 trees)
-   Map\_mammals\_median-LexpPD-p1.grd, .gri (gridded map of median expected loss in PD if all species present in a grid cell became extinct over 100 trees)
-   Map\_mammals\_LexpPD\_hotspots.grd (raster of loss-significant areas = the 2.5% cells with the highest expected loss in PD if all species present in a grid cell became extinct)
-   Map\_mammals\_prop-richness-threatened.grd, .gri (gridded map of the proportion of threatened species)
-   Map\_mammals\_prop-richness-rare.grd, .gri (gridded map of the proportion of rare species)
-   Map\_mammals\_mean-weighted-rarity.grd, .gri (gridded map of mean species-weighted rarity)
-   Map\_mammals\_prop-threatened-PD.grd, .gri (gridded map of the proportion of threatened PD)
-   Map\_mammals\_prop-rare-PD.grd, .gri (gridded map of the proportion of rare PD)
-   Map\_mammals\_mean-PWR.grd, .gri (gridded map of mean phylogenetic-weighted rarity)
-   Map\_mammals\_prop-richness-TOPHEDGE.grd, .gri (gridded map of the proportion of TOP 25% HEDGE species)
-   Map\_mammals\_prop-richness-TOPLEDGE.grd, .gri (gridded map of the proportion of TOP 25% LEDGE species)
-   tableS1.csv (Table S1 of the paper)
-   fig1.png (Figure 1 of the paper)
-   figS1.png (Figure S1 of the paper)
-   tableS2.csv (Table S2 of the paper)
-   tableS3.csv (Table S3 of the paper)
-   tableS4.csv (Table S4 of the paper)
-   fig2.png (Figure 2 of the paper)
-   fig3.png (Figure 3 of the paper)
-   figS2.png (Figure S2 of the paper)
-   fig4.png (Figure 4 of the paper)
-   figS3.png (Figure S3 of the paper)
-   tableS5.csv (Table S5 of the paper)
-   fig5.png (Figure 5 of the paper)
-   figS4.png (Figure S4 of the paper)
-   fig6.png (Figure 6 of the paper)

2. Requirements
---------------

Before starting the analyses, make sure to follow the three requirements described below.

### 2.1. clone repository and create additional folders

Clone this repository in your R working directory on your computer. This will create a folder "consmampd" in your R working directory with 4 folders (data, images, scripts, outputs). Then, in "consmampd", create one additional folder and name it "requireddata".

### 2.2. download data

Go to <https://datadryad.org/stash/dataset/doi:10.5061/dryad.bp26v20> and download data. Decompress the downloaded archive folder into the folder "requireddata". Repeat the operation for the .zip file "Current".

Download the supplementary files MAMMALS.shp, .shx, .sbx, .sbn, .prj, .dbf, .cpg and place them into the folder "requireddata". QUESTION: where to we put these big files (1.45 Go)? Dryad?

### 2.3. install packages

Install the following R packages:

``` r
install.packages(c("ape", "phylobase", "raster", "redlistr", "rgdal", "adiv", "rredlist", "picante", "phytools", "foreach", "doParallel", "Rarity", "dplyr", "ggplot2", "RColorBrewer", "tmap", "png", "grid"))
```

3. Analyses
-----------

Follow the three steps below to reproduce the analyses of the paper.

### 3.1. step 1: calculate species scores and identify species of interest

Run the script "calc\_species\_scores.R". This will allow you to:

-   calculate the HEDGE, LEDGE and ED scores for species (note that other scores are calculated but not used in the paper)
-   identify the TOP 25% HEDGE species, list the conservation measures they are benefitting and assign them an introduction status
-   identify the TOP 25% LEDGE species, list the conservation measures they are benefitting and assign them an introduction status.

### 3.2. step 2: calculate cells scores and identify areas of interest

Run the script "calc\_cells\_scores.R". This will allow you to:

-   calculate spatial scores for each cell of the reference gridded map
-   identify the 2.5% cells with the highest expected gain in PD if all species present in a grid cell were secured
-   identify the 2.5% cells with the highest expected loss in PD if all species present in a grid cell became extinct.

### 3.3. step 3: generate figures and tables

Run the script "figures\_tables.R". This will allow you to generate all the figures and tables of the paper. Note that, although this script is the backbone to generate the figures and tables of the paper, the final figures and tables appearing in the paper may slightly differ from those generated with this script because of customizations out of R.
