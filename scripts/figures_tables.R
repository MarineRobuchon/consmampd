###############################################################################################
# script to generate all figures and tables of the paper
###############################################################################################

### load packages and functions ###----
source("./scripts/new_functions.R")
library(RVAideMemoire)
library(ggplot2)
library(dplyr)
library(png)
library(grid)
library(raster)
library(RColorBrewer)
library(tmap)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


### Table S1 - ED, real and imputed status, GE, EDGE, LEDGE and HEDGE scores ----
mammals_scores <- read.csv2(paste0(getwd(), "/outputs/mammals25092019.csv"))[, -1]
mammals_traits <- read.csv2(paste0(getwd(), "/outputs/Trait_data_status_imputed.csv"))
mammals <- merge(mammals_scores, mammals_traits[, c("Binomial.1.2", "Order.1.2", "Family.1.2", "Imputed.Status")], by.x = "species_phylacine", by.y = "Binomial.1.2")
mammals$ED_rank <- rank(-mammals$ED, ties.method = "average")
mammals$GE <- status2pext(mammals$Imputed.Status)$GE
mammals$EDGE <- log(1+mammals$ED) + mammals$GE * log(2)
mammals$EDGE_rank <- rank(-mammals$EDGE, ties.method = "average")

tableS1 <- mammals[, c("species_phylacine", "Order.1.2.y", "Family.1.2.y", "IUCN.Status.1.2", "Imputed.Status", "GE", "ED", "ED_rank",
                       "EDGE", "EDGE_rank", "hedge", "hedge_rank", "ledge", "ledge_rank")]

colnames(tableS1) <- c("species", "order", "family", "IUCN_status_real", "IUCN_status_imputed", "GE", "ED","ED_rank",
                       "EDGE", "EDGE_rank", "HEDGE", "HEDGE_rank", "LEDGE", "LEDGE_rank")

write.csv2(tableS1, paste0(getwd(), "/outputs/tableS1.csv"))

### Figure 2 - Illustrating species in TOP 10 HEDGE (a) and TOP 10 LEDGE (b) ----
## TOP 10 HEDGE
top10hedge <- mammals[mammals$hedge_rank < 11 , c("species_phylacine", "hedge", "hedge_rank", "IUCN.Status.1.2")]
top10hedge <- top10hedge[order(top10hedge$hedge_rank),]
top10hedge$species <- gsub("_", " ", top10hedge$species_phylacine)
top10hedge$species[6] <- "Daubentonia\nmadagascariensis"
top10hedge$species[7] <- "Dicerorhinus\nsumatrensis"
top10hedge$species[8] <- "Solenodon\ncubanus"
top10hedge$species[9] <- "Solenodon\nparadoxus"
colnames(top10hedge)[4] <- "Red List status"
top10hedge$`Red List status` <- factor(top10hedge$`Red List status`)
levels(top10hedge$`Red List status`) <- c("EN", "CR")
top10hedge$species <- factor(top10hedge$species, levels = top10hedge$species)
top10hedge$species <- factor(top10hedge$species, levels = rev(levels(top10hedge$species)))
top10hedge$images <- list.files(path = paste0(getwd(), "/images/top10hedge"), full.names = TRUE)
top10hedge$reverse <- c(10:1)

plot <- ggplot(top10hedge, aes(x = species, y = hedge, fill = `Red List status`)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(values = c("orange", "red")) +
  scale_x_discrete(labels = rev(c("1.0", "2.0", "3.0", "4.5", "4.5", "6.0", "7.0", "8.5", "8.5", "10.0"))) +
  geom_text(aes(label = species), position = position_stack(vjust = 0.5), fontface = "italic", size = 3) +
  labs(title = "(a) TOP 10 HEDGE species", 
       x = "", y = "HEDGE (Ma)") +
  coord_flip() +
  theme_bw() +
  theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 9)) +
  scale_y_continuous(limits = c(0,85))

gr_hedge <- list()
for(i in 1:nrow(top10hedge)){
  img <-  readPNG(top10hedge$images[i])
  gr_hedge[[i]] <-   rasterGrob(img, interpolate=TRUE)
  plot <-  plot +
    annotation_custom(grob = gr_hedge[[i]], xmin = top10hedge$reverse[i]-.75, xmax = top10hedge$reverse[i]+.75, ymin = top10hedge$hedge[i]+0.5, ymax = top10hedge$hedge[i]+8)
}

fig2a <- plot
fig2a


## TOP 10 LEDGE
top10ledge <- mammals[mammals$ledge_rank < 11 , c("species_phylacine", "ledge", "ledge_rank", "IUCN.Status.1.2")]
top10ledge <- top10ledge[order(top10ledge$ledge_rank),]
top10ledge$species <- gsub("_", " ", top10ledge$species_phylacine)
colnames(top10ledge)[4] <- "Red List status"
top10ledge$`Red List status` <- factor(top10ledge$`Red List status`)
levels(top10ledge$`Red List status`) <- c("LC", "NT", "VU")
top10ledge$species <- factor(top10ledge$species, levels = top10ledge$species)
top10ledge$species <- factor(top10ledge$species, levels = rev(levels(top10ledge$species)))
top10ledge$images <- list.files(path = paste0(getwd(), "/images/top10ledge"), full.names = TRUE)
top10ledge$reverse <- c(10:1)

plot <- ggplot(top10ledge, aes(x = species, y = ledge, fill = `Red List status`)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_x_discrete(labels = rev(paste0(top10ledge$ledge_rank, ".0"))) +
  geom_text(aes(label = species), position = position_stack(vjust = 0.5), fontface = "italic", size = 3) +
  scale_fill_manual(values = c("green4", "greenyellow", "yellow")) +
  labs(title = "(b) TOP 10 LEDGE species", 
       x = "", y = "LEDGE (Ma)") +
  coord_flip() +
  theme_bw() +
  theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 9)) +
  scale_y_continuous(limits = c(0,85)) 

gr_ledge <- list()
for(i in 1:nrow(top10ledge)){
  img <-  readPNG(top10ledge$images[i])
  gr_ledge[[i]] <-   rasterGrob(img, interpolate=TRUE)
  plot <-  plot +
    annotation_custom(grob = gr_ledge[[i]], xmin = top10ledge$reverse[i]-.75, xmax = top10ledge$reverse[i]+.75, ymin = top10ledge$ledge[i]+0.5, ymax = top10ledge$ledge[i]+8)
}

fig2b <- plot
fig2b

## save Figure 2
ggsave(paste0(getwd(),"/outputs/fig2.png"), plot = multiplot(fig2a, fig2b, cols = 1), scale = 1, 
       width = 17, height = 25, units = "cm", dpi = 600, limitsize = TRUE)


### Figure S3 - Pairwise correlations between EDGE, HEDGE and LEDGE scores ----
png(filename = paste0(getwd(), "/outputs/figS3.png"), width = 17, height = 15, units = "cm", res = 600)
corPlot(df = tableS1[, c("EDGE", "HEDGE", "LEDGE")], method = "spearman", digits = 2, ties.method =  "average")
dev.off()

### Table S2 - Common species between TOP EDGE, TOP HEDGE, and TOP LEDGE ----
# comparing TOP HEDGE and TOP EDGE
tophedge <- head(mammals[order(mammals$hedge_rank) ,], n = round(nrow(mammals)*25/100))
topedge <- head(mammals[order(mammals$EDGE_rank) ,], n = round(nrow(mammals)*25/100))
species_hedge_edge <- intersect(tophedge$species_phylacine, topedge$species_phylacine)
length(species_hedge_edge) # 1194 common species
length(species_hedge_edge)/nrow(tophedge) *100 # 87 %
species_hedge_notedge <- setdiff(tophedge$species_phylacine, topedge$species_phylacine)
species_hedge_notedge # 175 species TOP HEDGE not in TOP EDGE
hedge_notedge <- tophedge[tophedge$species_phylacine %in% species_hedge_notedge,]

# comparing TOP LEDGE and TOP HEDGE
topledge <- head(mammals[order(mammals$ledge_rank) ,], n = round(nrow(mammals)*25/100))
species_ledge_hedge <- intersect(topledge$species_phylacine, tophedge$species_phylacine)
length(species_ledge_hedge) # 269 common species
length(species_ledge_hedge)/nrow(topledge) *100 # 20 %
ledge_hedge <- topledge[topledge$species_phylacine %in% species_ledge_hedge,]

# comparing TOP LEDGE and TOP EDGE
species_ledge_edge <- intersect(topledge$species_phylacine, topedge$species_phylacine)
length(species_ledge_edge) # 422 common species
length(species_ledge_edge)/nrow(topledge) *100 # 31 %
ledge_edge <- topledge[topledge$species_phylacine %in% species_ledge_edge,]

tableS2 <- data.frame(row.names = c("EDGE", "HEDGE", "LEDGE"))
tableS2$EDGE <- c(100,
                  length(species_hedge_edge)/nrow(tophedge)*100,
                  length(species_ledge_edge)/nrow(topledge)*100)
tableS2$HEDGE <- c(NA,
                   100,
                   length(species_ledge_hedge)/nrow(topledge)*100)
tableS2$LEDGE <- c(NA,
                   NA,
                   100)

write.csv2(tableS2, paste0(getwd(), "/outputs/tableS2.csv"))

### Table S4 - TOP 25% HEDGE: conservation and introduction ----
alldata_tophedge <- read.csv2(paste0(getwd(), "/outputs/mammals_tophedge_conservation_introduction.csv"))
tableS4 <- unique(alldata_tophedge[, c("species_phylacine", "hedge_rank", "hedge", "conservation_classification", "introduced")])
tableS4 <- tableS4[order(tableS4$hedge_rank),]
colnames(tableS4) <- c("species", "HEDGE_rank", "HEDGE", "conservation_classification", "introduced")
write.csv2(tableS4, paste0(getwd(), "/outputs/tableS4.csv"))

# species not protected
tophedge_notprotected <- tableS4[tableS4$conservation_classification=="None",]
nrow(tophedge_notprotected)/nrow(tophedge)*100

# species introduced
tophedge_introduced <- tableS4[tableS4$introduced=="yes",]
nlevels(factor(tophedge_introduced$species))/nrow(tophedge)*100

# species introduced AND on the list of Union concern
concern <- read.csv2(paste0(getwd(), "/data/animal_species_Union-concern.csv"))
intersect(gsub("_", " ", tophedge_introduced$species), concern$Scientific.name) # no tophedge species is on the list of Union concern

# species by realm
nrow(tophedge[tophedge$Marine=="1",])/nrow(mammals[mammals$Marine=="1",]) # 22% of marine mammals are TOP HEDGE
nrow(tophedge[tophedge$Terrestrial=="1",])/nrow(mammals[mammals$Terrestrial=="1",]) # 27% of terrestrial mammals are TOP HEDGE
nrow(tophedge[tophedge$Freshwater=="1",])/nrow(mammals[mammals$Freshwater=="1",]) # 41% of freshwater mammals are TOP HEDGE
nrow(tophedge[tophedge$Aerial=="1",])/nrow(mammals[mammals$Aerial=="1",]) # 18% of aerial mammals are TOP HEDGE

### Table S5 - TOP 25% LEDGE: conservation and introduction ----
alldata_topledge <- read.csv2(paste0(getwd(), "/outputs/mammals_topledge_conservation_introduction.csv"))
tableS5 <- unique(alldata_topledge[, c("species_phylacine", "ledge_rank", "ledge", "conservation_classification", "introduced")])
tableS5 <- tableS5[order(tableS5$ledge_rank),]
colnames(tableS5) <- c("species", "LEDGE_rank", "LEDGE", "conservation_classification", "introduced")
write.csv2(tableS5, paste0(getwd(), "/outputs/tableS5.csv"))

# species not protected
topledge_notprotected <- tableS5[tableS5$conservation_classification=="None",]
nrow(topledge_notprotected)/nrow(topledge)*100

# species introduced
topledge_introduced <- tableS5[tableS5$introduced=="yes",]
nlevels(factor(topledge_introduced$species))/nrow(topledge)*100

# species introduced that belong to both TOP LEDGE and TOP HEDGE
species_topledge_introduced <- unique(as.character(topledge_introduced$species)) # 34 species
species_tophedge_introduced <- unique(as.character(tophedge_introduced$species)) # 31 species
species_topledge_introduced[species_topledge_introduced %in% species_tophedge_introduced] # the 10 species species introduced that belong to both TOP LEDGE and TOP HEDGE

# species introduced AND on the list of Union concern
intersect(gsub("_", " ", topledge_introduced$species), concern$Scientific.name) # 2 top ledge species are on the list of Union concern

# species by realm
nrow(topledge[topledge$Marine=="1",])/nrow(mammals[mammals$Marine=="1",]) # 27% of marine mammals are TOP LEDGE
nrow(topledge[topledge$Terrestrial=="1",])/nrow(mammals[mammals$Terrestrial=="1",]) # 25% of terrestrial mammals are TOP LEDGE
nrow(topledge[topledge$Freshwater=="1",])/nrow(mammals[mammals$Freshwater=="1",]) # 34% of freshwater mammals are TOP LEDGE
nrow(topledge[topledge$Aerial=="1",])/nrow(mammals[mammals$Aerial=="1",]) # 25% of aerial mammals are TOP LEDGE



### Figure 5 ----
# Conservation measures for the TOP 25% HEDGE species (a) 
# and for the TOP 25% LEDGE species (b)
## fig5a
alldata_tophedge <- read.csv2(paste0(getwd(), "/outputs/mammals_tophedge_conservation_introduction.csv"))
gg_conservation_tophedge <- unique(alldata_tophedge[, c("species_phylacine", "conservation_classification", "introduced")])
colnames(gg_conservation_tophedge) <- c("species", "conservation", "introduced")
gg_conservation_tophedge$conservation <- as.factor(gg_conservation_tophedge$conservation)
# rename the levels
levels(gg_conservation_tophedge$conservation) <- c("Land/water protection", 
                                                   "Land/water management", 
                                                   "Species management", 
                                                   "Education & awareness", 
                                                   "Law & policy", 
                                                   "Livelihood, economic & other incentives", 
                                                   "None")
# re-order the levels so they fit the number of conservation measures
gg_conservation_tophedge$conservation <- factor(gg_conservation_tophedge$conservation, 
                                                levels = c ("Livelihood, economic & other incentives",
                                                            "Law & policy",
                                                            "Education & awareness",
                                                            "Species management",
                                                            "None",
                                                            "Land/water protection",
                                                            "Land/water management"))
# format the table to count the number of species by conservation and introduced
tab_conservation_tophedge <- gg_conservation_tophedge %>% count (conservation, introduced, 
                                                                 sort = TRUE)
# prepare the images
hedge05 <-  rasterGrob(readPNG(paste0(getwd(), "/images/top10hedge/hedge05.png")), interpolate = TRUE)
hedge103 <-  rasterGrob(readPNG(paste0(getwd(), "/images/introduced/hedge103.png")), interpolate = TRUE)
hedge02 <- rasterGrob(readPNG(paste0(getwd(), "/images/top10hedge/hedge02.png")), interpolate = TRUE)
# make the figure
fig5a <- ggplot(tab_conservation_tophedge, aes(x = conservation, y = n, fill = introduced)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#2166ac", "#8c510a")) +
  labs(title = "(a) Conservation measures for species in the TOP 25% HEDGE", 
       x = "Type of conservation measure", y = "Number of species in the TOP 25% HEDGE") +
  coord_flip(clip = "off") +
  theme_bw() +
  theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 9),
        axis.text.y = element_text(angle = 15)) +
  scale_y_continuous(limits = c(-550,1200), breaks = c(0, 250, 500, 750, 1000)) +
  annotation_custom(hedge05, xmin = 6.5, xmax = 7.5, ymin = -Inf, ymax = -340) +
  annotation_custom(hedge05, xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = -340) +
  annotation_custom(hedge05, xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = -340) +
  annotation_custom(hedge05, xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = -340) +
  annotate("text", x = 7, y = -180, label = "Varecia\nvariegata", fontface = "italic", size = 2, col = "#8c510a") +
  annotate("text", x = 6, y = -180, label = "Varecia\nvariegata", fontface = "italic", size = 2, col = "#8c510a") +
  annotate("text", x = 4, y = -180, label = "Varecia\nvariegata", fontface = "italic", size = 2, col = "#8c510a") +
  annotate("text", x = 3, y = -180, label = "Varecia\nvariegata", fontface = "italic", size = 2, col = "#8c510a") +
  annotation_custom(hedge103, xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = -340) +
  annotate("text", x = 5, y = -180, label = "Dasyprocta\nmexicana", fontface = "italic", size = 2, col = "#8c510a") +
  annotation_custom(hedge02, xmin = 4.5, xmax = 5.5, ymin = 460, ymax = 660) +
  annotate("text", x = 5, y = 780, label = "Mystacina\nrobusta", fontface = "italic", size = 2, col = "#2166ac")
fig5a

## fig5b
alldata_topledge <- read.csv2(paste0(getwd(), "/outputs/mammals_topledge_conservation_introduction.csv"))
gg_conservation_topledge <- unique(alldata_topledge[, c("species_phylacine", "conservation_classification", "introduced")])
colnames(gg_conservation_topledge) <- c("species", "conservation", "introduced")
gg_conservation_topledge$conservation <- as.factor(gg_conservation_topledge$conservation)
# rename the levels
levels(gg_conservation_topledge$conservation) <- c("Land/water protection", 
                                                   "Land/water management", 
                                                   "Species management", 
                                                   "Education & awareness", 
                                                   "Law & policy", 
                                                   "Livelihood, economic & other incentives", 
                                                   "None")
# re-order the levels so they fit the number of conservation measures
gg_conservation_topledge$conservation <- factor(gg_conservation_topledge$conservation, 
                                                levels = c ("Livelihood, economic & other incentives",
                                                            "Law & policy",
                                                            "Education & awareness",
                                                            "Species management",
                                                            "Land/water protection",
                                                            "Land/water management",
                                                            "None" ))
# format the table to count the number of species by conservation and introduced
tab_conservation_topledge <- gg_conservation_topledge %>% count (conservation, introduced, 
                                                                 sort = TRUE)
# prepare the images
ledge14 <-  rasterGrob(readPNG(paste0(getwd(), "/images/introduced/ledge14.png")), interpolate = TRUE)
ledge98 <-  rasterGrob(readPNG(paste0(getwd(), "/images/introduced/ledge98.png")), interpolate = TRUE)
ledge01 <- rasterGrob(readPNG(paste0(getwd(), "/images/top10ledge/ledge01.png")), interpolate = TRUE)
# make the figure
fig5b <- ggplot(tab_conservation_topledge, aes(x = conservation, y = n, fill = introduced)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#2166ac", "#8c510a")) +
  labs(title = "(b) Conservation measures for species in the TOP 25% LEDGE", 
       x = "Type of conservation measure", y = "Number of species in the TOP 25% LEDGE") +
  coord_flip(clip = "off") +
  theme_bw() +
  theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 9),
        axis.text.y = element_text(angle = 15)) +
  scale_y_continuous(limits = c(-550,1200), breaks = c(0, 250, 500, 750, 1000)) +
  annotation_custom(ledge14, xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = -340) +
  annotation_custom(ledge14, xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = -340) +
  annotation_custom(ledge14, xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = -340) +
  annotate("text", x = 5, y = -180, label = "Ornithorhynchus\nanatinus", fontface = "italic", size = 2, col = "#8c510a") +
  annotate("text", x = 6, y = -180, label = "Ornithorhynchus\nanatinus", fontface = "italic", size = 2, col = "#8c510a") +
  annotate("text", x = 4, y = -180, label = "Ornithorhynchus\nanatinus", fontface = "italic", size = 2, col = "#8c510a") +
  annotation_custom(ledge98, xmin = 6.5, xmax = 7.5, ymin = -Inf, ymax = -340) +
  annotate("text", x = 7, y = -180, label = "Myocastor\ncoypus", fontface = "italic", size = 2, col = "#8c510a") +
  annotation_custom(ledge01, xmin = 6.5, xmax = 7.5, ymin = 750, ymax = 1000) +
  annotate("text", x = 7, y = 1120, label = "Orycteropus\nafer", fontface = "italic", size = 2, col = "#2166ac")
  
fig5b

## save Figure 5
ggsave(paste0(getwd(),"/outputs/fig5.png"), 
       plot = multiplot(fig5a, fig5b, cols = 1), scale = 1, width = 17, height = 17, units = "cm",
       dpi = 600, limitsize = TRUE)

### Figures 3 and S2 - Spatial patterns and hotspots of ----
# the TOP 25% HEDGE species (a),
# the TOP 25% LEDGE species (b)

## 3a
richness_TOPHEDGE <- raster(paste0(getwd(), "/outputs/Map_mammals_richnessTOPHEDGE.grd"))
cuts <- quantile(richness_TOPHEDGE, c(0, 0.1, 0.6, 0.90, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- rgb(27, 158, 119, max = 255)

richness_TOPHEDGE.tm <- tm_shape(richness_TOPHEDGE) + 
  tm_raster (palette = colors, style = "fixed", title = "No. of species",
             breaks = cuts,
             labels = c("0", "1-3", "4-5","6-7", "8-9", "10-42")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(a) TOP 25% HEDGE species, 1369 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))

# load and transform the coastline spatial file
coastline <- shapefile(paste0(getwd(),"/data/ne_50m_coastline.shp"))
coastline.tm <- tm_shape(coastline, projection = crs(richness_TOPHEDGE)) + tm_lines(lwd = 0.5)
  
richness_TOPHEDGE.tm + coastline.tm 
fig3a <- richness_TOPHEDGE.tm + coastline.tm

## S2aprop
prop_richness_TOPHEDGE <- raster(paste0(getwd(), "/outputs/Map_mammals_prop-richness-TOPHEDGE.grd"))
cuts <- quantile(prop_richness_TOPHEDGE, c(0, 0.1, 0.6, 0.90, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- rgb(27, 158, 119, max = 255)

prop_richness_TOPHEDGE.tm <- tm_shape(prop_richness_TOPHEDGE) + 
  tm_raster (palette = colors, style = "fixed", title = "No. of TOP 25% HEDGE species/\nNo. of species (%)",
             breaks = cuts,
             labels = c("0-3.2", "3.3-14.9", "15.0-20.1","20.2-24.9", "25.0-29.3", "29.4-66.7")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(a) Proportion of TOP 25% HEDGE species, 1369 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
prop_richness_TOPHEDGE.tm + coastline.tm 
figS2aprop <- prop_richness_TOPHEDGE.tm + coastline.tm

## 3b
richness_TOPLEDGE <- raster(paste0(getwd(), "/outputs/Map_mammals_richnessTOPLEDGE.grd"))
cuts <- quantile(richness_TOPLEDGE, c(0, 0.1, 0.5, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- rgb(217, 95, 2, max = 255)

richness_TOPLEDGE.tm <- tm_shape(richness_TOPLEDGE) + 
  tm_raster (palette = colors, style = "fixed", title = "No. of species",
             breaks = cuts,
             labels = c("0-5", "6-7", "8-25", "26-40", "41-53", "54-98")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(b) TOP 25% LEDGE species, 1369 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
richness_TOPLEDGE.tm + coastline.tm 
fig3b <- richness_TOPLEDGE.tm + coastline.tm

## S2bprop
prop_richness_TOPLEDGE <- raster(paste0(getwd(), "/outputs/Map_mammals_prop-richness-TOPLEDGE.grd"))
cuts <- quantile(prop_richness_TOPLEDGE, c(0, 0.1, 0.5, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- rgb(217, 95, 2, max = 255)

prop_richness_TOPLEDGE.tm <- tm_shape(prop_richness_TOPLEDGE) + 
  tm_raster (palette = colors, style = "fixed", title = "No. of TOP 25% LEDGE species/\nNo. of species (%)",
             breaks = cuts,
             labels = c("0-27.7", "27.8-34.5", "34.6-43.8", "43.9-47.3", "47.4-49.9", "50.0-100.0")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(b) Proportion of TOP 25% LEDGE species, 1369 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
prop_richness_TOPLEDGE.tm + coastline.tm 
figS2bprop <- prop_richness_TOPLEDGE.tm + coastline.tm

## save Figure 3 
png(filename = paste0(getwd(), "/outputs/fig3.png"),
    width = 13, height = 8, units = "cm", pointsize = 8, res = 600)
multiplot(fig3a, fig3b, cols = 1)
dev.off()

## save Figure S2
png(filename = paste0(getwd(), "/outputs/figS1.png"),
    width = 13, height = 8, units = "cm", pointsize = 8, res = 600)
multiplot(figS2aprop, figS2bprop, cols = 1)
dev.off()

### Figure 4 - Spatial patterns and hotspots of ----
# expected gain in phylogenetic diversity if all species present in the cell are saved from extinction (p = 0) (a) 
# and expected loss in phylogenetic diversity if all species present in the cell become extinct (p = 1) (b)

## 4a 
median_GexpPD_p0 <- raster(paste0(getwd(), "/outputs/Map_mammals_median-GexpPD-p0.grd"))
cuts <- quantile(median_GexpPD_p0, c(0, 0.1, 0.5, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- rgb(27, 158, 119, max = 255)

median_GexpPD_p0.tm <- tm_shape(median_GexpPD_p0) + 
  tm_raster (palette = colors, style = "fixed", title = "Evolutionary history (Ma)",
             breaks = cuts, labels = c("0-0.4", "0.5-5.4", "5.5-7.9", "8.0-12.0", "12.1-17.2", "17.3-78.4")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(a) Gain in expected phylogenetic diversity if all species present in the cell are saved from extinction",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
median_GexpPD_p0.tm + coastline.tm 
fig4a <- median_GexpPD_p0.tm + coastline.tm 

## 4b
median_LexpPD_p1 <- raster(paste0(getwd(), "/outputs/Map_mammals_median-LexpPD-p1.grd"))
cuts <- quantile(abs(median_LexpPD_p1), c(0, 0.1, 0.5, 0.9, 0.95, 0.975, 1)) # abs because values are negative (it is a loss)
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- rgb(217, 95, 2, max = 255)

median_LexpPD_p1.tm <- tm_shape(abs(median_LexpPD_p1)) + 
  tm_raster (palette = colors, style = "fixed", title = "Evolutionary history (Ma)",
             breaks = cuts, labels = c("1.2-87.8", "87.9-105.6", "105.7-310.7", "310.8-532.9", "533.0-696.2", "696.3-1172.0")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(b) Loss in expected phylogenetic diversity if all species present in the cell become extinct",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
median_LexpPD_p1.tm + coastline.tm 
fig4b <- median_LexpPD_p1.tm + coastline.tm 

## Save Figure 4
png(filename = paste0(getwd(), "/outputs/fig4.png"),
    width = 13, height = 8, units = "cm", pointsize = 8, res = 600)
multiplot(fig4a, fig4b, cols = 1)
dev.off()

### Figure S4 - Pairwise correlations between the 12 spatial scores ----
SR <- raster(paste0(getwd(), "/outputs/Map_mammals_richness"))
TSR <- raster(paste0(getwd(), "/outputs/Map_mammals_threatened-richness"))
RSR <- raster(paste0(getwd(), "/outputs/Map_mammals_rare-richness"))
SWR <- raster(paste0(getwd(), "/outputs/Map_mammals_weighted-rarity"))
PD <- raster(paste0(getwd(), "/outputs/Map_mammals_median-PD"))
TPD <- raster(paste0(getwd(), "/outputs/Map_mammals_median-threatened-PD"))
RPD <- raster(paste0(getwd(), "/outputs/Map_mammals_median-rare-PD"))
PWR <- raster(paste0(getwd(), "/outputs/Map_mammals_median-PWR"))
HEDGE <- raster(paste0(getwd(), "/outputs/Map_mammals_richnessTOPHEDGE"))
LEDGE <- raster(paste0(getwd(), "/outputs/Map_mammals_richnessTOPLEDGE"))
GexpPD <- raster(paste0(getwd(), "/outputs/Map_mammals_median-GexpPD-p0"))
LexpPD <- raster(paste0(getwd(), "/outputs/Map_mammals_median-LexpPD-p1"))

table_spatial_scores <- data.frame(SR = getValues(SR), TSR = getValues(TSR), RSR = getValues(RSR), SWR = getValues(SWR),
                                   PD = getValues(PD), TPD = getValues(TPD), RPD = getValues(RPD), PWR = getValues(PWR),
                                   HEDGE = getValues(HEDGE), LEDGE = getValues(abs(LEDGE)), GexpPD = getValues(GexpPD), LexpPD = getValues(abs(LexpPD)))
head(table_spatial_scores)

png(filename = paste0(getwd(), "/outputs/figS4.png"), width = 17, height = 17, units = "cm", res = 600, pointsize = 8)
corPlot(df = table_spatial_scores, method = "spearman", digits = 2, ties.method =  "average")
dev.off()

### Table S3 - Common zones between hotspots for the 12 spatial scores ----
SR_hotspots <- SR
SR_hotspots[SR_hotspots < quantile(SR, 0.975)] <- 0
SR_hotspots[SR_hotspots >= quantile(SR, 0.975)] <- 1

TSR_hotspots <- TSR
TSR_hotspots[TSR_hotspots < quantile(TSR, 0.975)] <- 0
TSR_hotspots[TSR_hotspots >= quantile(TSR, 0.975)] <- 1

RSR_hotspots <- RSR
RSR_hotspots[RSR_hotspots < quantile(RSR, 0.975)] <- 0
RSR_hotspots[RSR_hotspots >= quantile(RSR, 0.975)] <- 1

SWR_hotspots <- SWR
SWR_hotspots[SWR_hotspots < quantile(SWR, 0.975)] <- 0
SWR_hotspots[SWR_hotspots >= quantile(SWR, 0.975)] <- 1

PD_hotspots <- PD
PD_hotspots[PD_hotspots < quantile(PD, 0.975)] <- 0
PD_hotspots[PD_hotspots >= quantile(PD, 0.975)] <- 1

TPD_hotspots <- TPD
TPD_hotspots[TPD_hotspots < quantile(TPD, 0.975)] <- 0
TPD_hotspots[TPD_hotspots >= quantile(TPD, 0.975)] <- 1

RPD_hotspots <- RPD
RPD_hotspots[RPD_hotspots < quantile(RPD, 0.975)] <- 0
RPD_hotspots[RPD_hotspots >= quantile(RPD, 0.975)] <- 1

PWR_hotspots <- PWR
PWR_hotspots[PWR_hotspots < quantile(PWR, 0.975)] <- 0
PWR_hotspots[PWR_hotspots >= quantile(PWR, 0.975)] <- 1

HEDGE_hotspots <- HEDGE
HEDGE_hotspots[HEDGE_hotspots < quantile(HEDGE, 0.975)] <- 0
HEDGE_hotspots[HEDGE_hotspots >= quantile(HEDGE, 0.975)] <- 1

LEDGE_hotspots <- LEDGE
LEDGE_hotspots[LEDGE_hotspots < quantile(LEDGE, 0.975)] <- 0
LEDGE_hotspots[LEDGE_hotspots >= quantile(LEDGE, 0.975)] <- 1

GexpPD_hotspots <- GexpPD
GexpPD_hotspots[GexpPD_hotspots < quantile(GexpPD, 0.975)] <- 0
GexpPD_hotspots[GexpPD_hotspots >= quantile(GexpPD, 0.975)] <- 1

LexpPD_hotspots <- abs(LexpPD)
LexpPD_hotspots[LexpPD_hotspots < quantile(abs(LexpPD), 0.975)] <- 0
LexpPD_hotspots[LexpPD_hotspots >= quantile(abs(LexpPD), 0.975)] <- 1

table_spatial_hotspots <- data.frame(SR = getValues(SR_hotspots), TSR = getValues(TSR_hotspots), 
                                     RSR = getValues(RSR_hotspots), SWR = getValues(SWR_hotspots),
                                     PD = getValues(PD_hotspots), TPD = getValues(TPD_hotspots), 
                                     RPD = getValues(RPD_hotspots), PWR = getValues(PWR_hotspots),
                                     HEDGE = getValues(HEDGE_hotspots), LEDGE = getValues(LEDGE_hotspots), 
                                     GexpPD = getValues(GexpPD_hotspots), LexpPD = getValues(LexpPD_hotspots))

tableS3 <- data.frame(row.names = c("SR", "TSR", "RSR", "SWR",
                                    "PD", "TPD", "RPD", "PWR",
                                    "HEDGE", "LEDGE", "GexpPD", "LexpPD"))

tableS3$SR <- c(100, 
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$SR == 1 & table_spatial_hotspots$TSR == 1) ,])/max(colSums(table_spatial_hotspots)[1], colSums(table_spatial_hotspots)[2])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$SR == 1 & table_spatial_hotspots$RSR == 1) ,])/max(colSums(table_spatial_hotspots)[1], colSums(table_spatial_hotspots)[3])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$SR == 1 & table_spatial_hotspots$SWR == 1) ,])/max(colSums(table_spatial_hotspots)[1], colSums(table_spatial_hotspots)[4])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$SR == 1 & table_spatial_hotspots$PD == 1) ,])/max(colSums(table_spatial_hotspots)[1], colSums(table_spatial_hotspots)[5])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$SR == 1 & table_spatial_hotspots$TPD == 1) ,])/max(colSums(table_spatial_hotspots)[1], colSums(table_spatial_hotspots)[6])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$SR == 1 & table_spatial_hotspots$RPD == 1) ,])/max(colSums(table_spatial_hotspots)[1], colSums(table_spatial_hotspots)[7])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$SR == 1 & table_spatial_hotspots$PWR == 1) ,])/max(colSums(table_spatial_hotspots)[1], colSums(table_spatial_hotspots)[8])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$SR == 1 & table_spatial_hotspots$HEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[1], colSums(table_spatial_hotspots)[9])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$SR == 1 & table_spatial_hotspots$LEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[1], colSums(table_spatial_hotspots)[10])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$SR == 1 & table_spatial_hotspots$GexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[1], colSums(table_spatial_hotspots)[11])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$SR == 1 & table_spatial_hotspots$LexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[1], colSums(table_spatial_hotspots)[12])*100)

tableS3$TSR <- c(NA, 
                 100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TSR == 1 & table_spatial_hotspots$RSR == 1) ,])/max(colSums(table_spatial_hotspots)[2], colSums(table_spatial_hotspots)[3])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TSR == 1 & table_spatial_hotspots$SWR == 1) ,])/max(colSums(table_spatial_hotspots)[2], colSums(table_spatial_hotspots)[4])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TSR == 1 & table_spatial_hotspots$PD == 1) ,])/max(colSums(table_spatial_hotspots)[2], colSums(table_spatial_hotspots)[5])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TSR == 1 & table_spatial_hotspots$TPD == 1) ,])/max(colSums(table_spatial_hotspots)[2], colSums(table_spatial_hotspots)[6])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TSR == 1 & table_spatial_hotspots$RPD == 1) ,])/max(colSums(table_spatial_hotspots)[2], colSums(table_spatial_hotspots)[7])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TSR == 1 & table_spatial_hotspots$PWR == 1) ,])/max(colSums(table_spatial_hotspots)[2], colSums(table_spatial_hotspots)[8])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TSR == 1 & table_spatial_hotspots$HEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[2], colSums(table_spatial_hotspots)[9])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TSR == 1 & table_spatial_hotspots$LEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[2], colSums(table_spatial_hotspots)[10])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TSR == 1 & table_spatial_hotspots$GexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[2], colSums(table_spatial_hotspots)[11])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TSR == 1 & table_spatial_hotspots$LexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[2], colSums(table_spatial_hotspots)[12])*100)

tableS3$RSR <- c(NA, 
                 NA,
                 100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RSR == 1 & table_spatial_hotspots$SWR == 1) ,])/max(colSums(table_spatial_hotspots)[3], colSums(table_spatial_hotspots)[4])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RSR == 1 & table_spatial_hotspots$PD == 1) ,])/max(colSums(table_spatial_hotspots)[3], colSums(table_spatial_hotspots)[5])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RSR == 1 & table_spatial_hotspots$TPD == 1) ,])/max(colSums(table_spatial_hotspots)[3], colSums(table_spatial_hotspots)[6])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RSR == 1 & table_spatial_hotspots$RPD == 1) ,])/max(colSums(table_spatial_hotspots)[3], colSums(table_spatial_hotspots)[7])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RSR == 1 & table_spatial_hotspots$PWR == 1) ,])/max(colSums(table_spatial_hotspots)[3], colSums(table_spatial_hotspots)[8])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RSR == 1 & table_spatial_hotspots$HEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[3], colSums(table_spatial_hotspots)[9])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RSR == 1 & table_spatial_hotspots$LEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[3], colSums(table_spatial_hotspots)[10])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RSR == 1 & table_spatial_hotspots$GexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[3], colSums(table_spatial_hotspots)[11])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RSR == 1 & table_spatial_hotspots$LexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[3], colSums(table_spatial_hotspots)[12])*100)

tableS3$SWR <- c(NA, 
                 NA,
                 NA,
                 100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$SWR == 1 & table_spatial_hotspots$PD == 1) ,])/max(colSums(table_spatial_hotspots)[4], colSums(table_spatial_hotspots)[5])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$SWR == 1 & table_spatial_hotspots$TPD == 1) ,])/max(colSums(table_spatial_hotspots)[4], colSums(table_spatial_hotspots)[6])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$SWR == 1 & table_spatial_hotspots$RPD == 1) ,])/max(colSums(table_spatial_hotspots)[4], colSums(table_spatial_hotspots)[7])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$SWR == 1 & table_spatial_hotspots$PWR == 1) ,])/max(colSums(table_spatial_hotspots)[4], colSums(table_spatial_hotspots)[8])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$SWR == 1 & table_spatial_hotspots$HEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[4], colSums(table_spatial_hotspots)[9])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$SWR == 1 & table_spatial_hotspots$LEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[4], colSums(table_spatial_hotspots)[10])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$SWR == 1 & table_spatial_hotspots$GexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[4], colSums(table_spatial_hotspots)[11])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$SWR == 1 & table_spatial_hotspots$LexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[4], colSums(table_spatial_hotspots)[12])*100)

tableS3$PD <- c(NA, 
                NA,
                NA,
                NA,
                100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$PD == 1 & table_spatial_hotspots$TPD == 1) ,])/max(colSums(table_spatial_hotspots)[5], colSums(table_spatial_hotspots)[6])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$PD == 1 & table_spatial_hotspots$RPD == 1) ,])/max(colSums(table_spatial_hotspots)[5], colSums(table_spatial_hotspots)[7])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$PD == 1 & table_spatial_hotspots$PWR == 1) ,])/max(colSums(table_spatial_hotspots)[5], colSums(table_spatial_hotspots)[8])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$PD == 1 & table_spatial_hotspots$HEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[5], colSums(table_spatial_hotspots)[9])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$PD == 1 & table_spatial_hotspots$LEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[5], colSums(table_spatial_hotspots)[10])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$PD == 1 & table_spatial_hotspots$GexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[5], colSums(table_spatial_hotspots)[11])*100,
                nrow(table_spatial_hotspots[which(table_spatial_hotspots$PD == 1 & table_spatial_hotspots$LexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[5], colSums(table_spatial_hotspots)[12])*100)

tableS3$TPD <- c(NA, 
                 NA,
                 NA,
                 NA,
                 NA,
                 100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TPD == 1 & table_spatial_hotspots$RPD == 1) ,])/max(colSums(table_spatial_hotspots)[6], colSums(table_spatial_hotspots)[7])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TPD == 1 & table_spatial_hotspots$PWR == 1) ,])/max(colSums(table_spatial_hotspots)[6], colSums(table_spatial_hotspots)[8])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TPD == 1 & table_spatial_hotspots$HEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[6], colSums(table_spatial_hotspots)[9])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TPD == 1 & table_spatial_hotspots$LEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[6], colSums(table_spatial_hotspots)[10])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TPD == 1 & table_spatial_hotspots$GexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[6], colSums(table_spatial_hotspots)[11])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$TPD == 1 & table_spatial_hotspots$LexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[6], colSums(table_spatial_hotspots)[12])*100)

tableS3$RPD <- c(NA, 
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RPD == 1 & table_spatial_hotspots$PWR == 1) ,])/max(colSums(table_spatial_hotspots)[7], colSums(table_spatial_hotspots)[8])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RPD == 1 & table_spatial_hotspots$HEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[7], colSums(table_spatial_hotspots)[9])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RPD == 1 & table_spatial_hotspots$LEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[7], colSums(table_spatial_hotspots)[10])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RPD == 1 & table_spatial_hotspots$GexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[7], colSums(table_spatial_hotspots)[11])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$RPD == 1 & table_spatial_hotspots$LexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[7], colSums(table_spatial_hotspots)[12])*100)

tableS3$PWR <- c(NA, 
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 NA,
                 100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$PWR == 1 & table_spatial_hotspots$HEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[8], colSums(table_spatial_hotspots)[9])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$PWR == 1 & table_spatial_hotspots$LEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[8], colSums(table_spatial_hotspots)[10])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$PWR == 1 & table_spatial_hotspots$GexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[8], colSums(table_spatial_hotspots)[11])*100,
                 nrow(table_spatial_hotspots[which(table_spatial_hotspots$PWR == 1 & table_spatial_hotspots$LexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[8], colSums(table_spatial_hotspots)[12])*100)

tableS3$HEDGE <- c(NA, 
                   NA,
                   NA,
                   NA,
                   NA,
                   NA,
                   NA,
                   NA,
                   100,
                   nrow(table_spatial_hotspots[which(table_spatial_hotspots$HEDGE == 1 & table_spatial_hotspots$LEDGE == 1) ,])/max(colSums(table_spatial_hotspots)[9], colSums(table_spatial_hotspots)[10])*100,
                   nrow(table_spatial_hotspots[which(table_spatial_hotspots$HEDGE == 1 & table_spatial_hotspots$GexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[9], colSums(table_spatial_hotspots)[11])*100,
                   nrow(table_spatial_hotspots[which(table_spatial_hotspots$HEDGE == 1 & table_spatial_hotspots$LexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[9], colSums(table_spatial_hotspots)[12])*100)

tableS3$LEDGE <- c(NA, 
                   NA,
                   NA,
                   NA,
                   NA,
                   NA,
                   NA,
                   NA,
                   NA,
                   100,
                   nrow(table_spatial_hotspots[which(table_spatial_hotspots$LEDGE == 1 & table_spatial_hotspots$GexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[10], colSums(table_spatial_hotspots)[11])*100,
                   nrow(table_spatial_hotspots[which(table_spatial_hotspots$LEDGE == 1 & table_spatial_hotspots$LexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[10], colSums(table_spatial_hotspots)[12])*100)

tableS3$GexpPD <- c(NA, 
                    NA,
                    NA,
                    NA,
                    NA,
                    NA,
                    NA,
                    NA,
                    NA,
                    NA,
                    100,
                    nrow(table_spatial_hotspots[which(table_spatial_hotspots$GexpPD == 1 & table_spatial_hotspots$LexpPD == 1) ,])/max(colSums(table_spatial_hotspots)[11], colSums(table_spatial_hotspots)[12])*100)

tableS3$LexpPD <- c(NA, 
                    NA,
                    NA,
                    NA,
                    NA,
                    NA,
                    NA,
                    NA,
                    NA,
                    NA,
                    NA,
                    100)

write.csv2(tableS3, paste0(getwd(), "/outputs/tableS3.csv"))

### Figure S5 and S6 - Spatial patterns and hotspots of ---- FIX§§§
# species richness (a), 
# threatened species richness (b), 
# rare species richness (c), 
# species-weighted rarity (d),
# phylogenetic diversity (e), 
# threatened phylogenetic diversity (f), 
# rare phylogenetic diversity (g),
# and phylogenetic-weighted rarity (h) 

## S5a
richness <- raster(paste0(getwd(), "/outputs/Map_mammals_richness.grd"))
cuts <- quantile(richness, c(0, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975, 1))
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"
richness.tm <- tm_shape(richness) + 
  tm_raster (palette = colors, style = "fixed", title = "No. of species",
             breaks = cuts,
             labels = c("1-20", "21-22", "23-32", "33-72", "73-119", "120-148", "149-242")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(a) Species richness, 5477 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
richness.tm + coastline.tm 
figS5a <- richness.tm + coastline.tm 

## S5b
richness_threatened <- raster(paste0(getwd(), "/outputs/Map_mammals_threatened-richness.grd"))
cuts <- quantile(richness_threatened, c(0, 0.1, 0.75, 0.90, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"

richness_threatened.tm <- tm_shape(richness_threatened) + 
  tm_raster (palette = colors, style = "fixed", title = "No. of species",
             breaks = cuts,
             labels = c("0", "1-3", "4-5", "6-7", "8-10", "11-44")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(b) Threatened species richness, 1193 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
richness_threatened.tm + coastline.tm 
figS5b <- richness_threatened.tm + coastline.tm

## S6bprop
prop_richness_threatened <- raster(paste0(getwd(), "/outputs/Map_mammals_prop-richness-threatened.grd"))
cuts <- quantile(prop_richness_threatened, c(0, 0.1, 0.75, 0.90, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"

prop_richness_threatened.tm <- tm_shape(prop_richness_threatened) + 
  tm_raster (palette = colors, style = "fixed", title = "No. of threatened species/\nNo. of species (%)",
             breaks = cuts,
             labels = c("0-3", "4-17", "18-20", "21-27", "28-32", "33-100")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(b) Proportion of threatened species, 1193 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
prop_richness_threatened.tm + coastline.tm 
figS6b <- prop_richness_threatened.tm + coastline.tm

## S5c
richness_rare <- raster(paste0(getwd(), "/outputs/Map_mammals_rare-richness.grd"))
cuts <- quantile(richness_rare, c(0, 0.85, 0.90, 0.95, 0.975, 1))
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"

richness_rare.tm <- tm_shape(richness_rare) + 
  tm_raster (palette = colors, style = "fixed", title = "No. of species",
             breaks = cuts, labels = c("0", "1", "2-3", "4-7", "8-57")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(c) Rare species richness, 2736 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
richness_rare.tm + coastline.tm 
figS5c <- richness_rare.tm + coastline.tm

## S6cprop
prop_richness_rare <- raster(paste0(getwd(), "/outputs/Map_mammals_prop-richness-rare.grd"))
cuts <- quantile(prop_richness_rare, c(0, 0.85, 0.90, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"

prop_richness_rare.tm <- tm_shape(prop_richness_rare) + 
  tm_raster (palette = colors, style = "fixed", title = "No. of rare species/\nNo. of species (%)",
             breaks = cuts, labels = c("< 1", "1-2", "3-5", "6-9", "10-81")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(c) Proportion of rare species, 2736 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
prop_richness_rare.tm + coastline.tm 
figS6c <- prop_richness_rare.tm + coastline.tm

## S5d
weighted_rarity <- raster(paste0(getwd(), "/outputs/Map_mammals_weighted-rarity.grd"))
cuts <- quantile(weighted_rarity, c(0, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de3d26"

weighted_rarity.tm <- tm_shape(weighted_rarity) + 
  tm_raster (palette = colors, style = "fixed", title = "Score",
             breaks = cuts, labels = c("0-0.2", "0.3-0.5", "0.6-0.7", "0.8-7.3")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(d) Species-weighted rarity, 5477 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
weighted_rarity.tm + coastline.tm
figS5d <- weighted_rarity.tm + coastline.tm

## S6dprop
prop_weighted_rarity <- raster(paste0(getwd(), "/outputs/Map_mammals_mean-weighted-rarity.grd"))
cuts <- quantile(prop_weighted_rarity, c(0, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de3d26"

prop_weighted_rarity.tm <- tm_shape(prop_weighted_rarity) + 
  tm_raster (palette = colors, style = "fixed", title = "Species-weighted rarity/\nSpecies richness",
             breaks = cuts, labels = c("0-0.003", "0.004-0.005", "0.006-0.008", "0.009-0.1")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(d) Mean species-weighted rarity, 5477 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
prop_weighted_rarity.tm + coastline.tm 
figS6d <- prop_weighted_rarity.tm + coastline.tm

## S5e
median_PD <- raster(paste0(getwd(), "/outputs/Map_mammals_median-PD.grd"))
cuts <- quantile(median_PD, c(0, 0.5, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"

median_PD.tm <- tm_shape(median_PD) + 
  tm_raster (palette = colors, style = "fixed", title = "Evolutionary history (Ma)",
             breaks = cuts, labels = c("217-503", "504-2193", "2194-3082", "3083-3485", "3486-4689")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(e) Phylogenetic diversity, 5477 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
median_PD.tm + coastline.tm 
figS5e <- median_PD.tm + coastline.tm

## S5f
median_threatened_PD <- raster(paste0(getwd(), "/outputs/Map_mammals_median-threatened-PD.grd"))
cuts <- quantile(median_threatened_PD, c(0, 0.1, 0.5, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"

median_threatened_PD.tm <- tm_shape(median_threatened_PD) + 
  tm_raster (palette = colors, style = "fixed", title = "Evolutionary history (Ma)",
             breaks = cuts, labels = c("0-217", "218-259", "260-455", "455-647", "648-820", "821-1793")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(f) Threatened phylogenetic diversity, 1193 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
median_threatened_PD.tm + coastline.tm 
figS5f <- median_threatened_PD.tm + coastline.tm

## S6fprop
prop_median_threatened_PD <- raster(paste0(getwd(), "/outputs/Map_mammals_prop-threatened-PD.grd"))
cuts <- quantile(prop_median_threatened_PD, c(0, 0.1, 0.5, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"

prop_median_threatened_PD.tm <- tm_shape(prop_median_threatened_PD) + 
  tm_raster (palette = colors, style = "fixed", title = "Threatened phylogenetic diversity/\nPhylogenetic diversity (%)",
             breaks = cuts, labels = c("0-15", "16-57", "58-64", "65-68", "69-72", "73-100")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(f) Proportion of threatened phylogenetic diversity, 1193 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
prop_median_threatened_PD.tm + coastline.tm 
figS6f<- prop_median_threatened_PD.tm + coastline.tm

## S5g
median_rare_PD <- raster(paste0(getwd(), "/outputs/Map_mammals_median-rare-PD.grd"))
cuts <- quantile(median_rare_PD, c(0, 0.85, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"

median_rare_PD.tm <- tm_shape(median_rare_PD) + 
  tm_raster (palette = colors, style = "fixed", title = "Evolutionary history (Ma)",
             breaks = cuts, labels = c("0-217", "218-245", "246-413", "414-586", "587-1523")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(g) Rare phylogenetic diversity, 2736 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
median_rare_PD.tm + coastline.tm 
figS5g <- median_rare_PD.tm + coastline.tm 

## S6gprop
prop_median_rare_PD <- raster(paste0(getwd(), "/outputs/Map_mammals_prop-rare-PD.grd"))
cuts <- quantile(prop_median_rare_PD, c(0, 0.85, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"

prop_median_rare_PD.tm <- tm_shape(prop_median_rare_PD) + 
  tm_raster (palette = colors, style = "fixed", title = "Rare phylogenetic diversity/\nPhylogenetic diversity (%)",
             breaks = cuts, labels = c("0-7.7", "7.8-13.1", "13.2-21.0", "21.1-28.1", "28.2-75.0")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(g) Proportion of rare phylogenetic diversity, 2736 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
prop_median_rare_PD.tm + coastline.tm 
figS6g <- prop_median_rare_PD.tm + coastline.tm 

## S5h
median_PWR <- raster(paste0(getwd(), "/outputs/Map_mammals_median-PWR.grd"))
cuts <- quantile(median_PWR, c(0, 0.5, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"

median_PWR.tm <- tm_shape(median_PWR) + 
  tm_raster (palette = colors, style = "fixed", title = "Evolutionary history (Ma)",
             breaks = cuts, labels = c("0-0.017", "0.018-0.200", "0.200-3.560","3.561-5.170", "5.171-46.656")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(h) Phylogenetic-weighted rarity, 5477 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
median_PWR.tm + coastline.tm
figS5h <- median_PWR.tm + coastline.tm

## S6hprop
prop_median_PWR <- raster(paste0(getwd(), "/outputs/Map_mammals_mean-PWR.grd"))
cuts <- quantile(prop_median_PWR, c(0, 0.9, 0.95, 0.975, 1))
cuts
colors <- brewer.pal(9, "Greys")
colors [9] <- "#de2d26"

prop_median_PWR.tm <- tm_shape(prop_median_PWR) + 
  tm_raster (palette = colors, style = "fixed", title = "Phylogenetic-weighted rarity/\nPhylogenetic diversity",
             breaks = cuts, labels = c("0-0.08", "0.09-0.13", "0.14-0.18","0.19-2.36")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right",
            main.title = "(h) Mean phylogenetic-weighted rarity, 5477 species",
            main.title.size = 1, outer.margins = c(0, 0, 0, 0))
prop_median_PWR.tm + coastline.tm
figS6h <- prop_median_PWR.tm + coastline.tm

## save Figure S5
png(filename = paste0(getwd(), "/outputs/figS5.png"),
    width = 17, height = 12, units = "cm", pointsize = 8, res = 600)
multiplot(figS5a, figS5b, figS5c, figS5d, figS5e, figS5f, figS5g, figS5h, cols = 2)
dev.off()

## save Figure S6
png(filename = paste0(getwd(), "/outputs/figS6.png"),
    width = 17, height = 12, units = "cm", pointsize = 8, res = 600)
multiplot(figS5a, figS6b, figS6c, figS6d, figS5e, figS6f, figS6g, figS6h, cols = 2)
dev.off()


### Figure 6 - Intersection between priority areas and and the percentage of area protected by the current network of protected areas ----
median_GexpPD_p0 <- raster(paste0(getwd(), "/outputs/Map_mammals_median-GexpPD-p0.grd"))
GexpPD_hotspots <- median_GexpPD_p0
GexpPD_hotspots[GexpPD_hotspots < quantile(median_GexpPD_p0, 0.975)] <- NA # NA values for cells which are not hotspots
GexpPD_hotspots[GexpPD_hotspots >= quantile(median_GexpPD_p0, 0.975)] <- 1 # assign value 1 to priority areas

ppa <-  raster(paste0(getwd(), "/data/percentage_pa_grid_raster.tif")) # raster with percentage of protected area in each grid cell

GexpPD_hotspots.tm <- tm_shape(GexpPD_hotspots) +
  tm_raster(palette = "Dark2", alpha = 0.9, labels = "", title = "priority areas") +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right", outer.margins = c(0, 0, 0, 0))

ppa.tm <- tm_shape(ppa) + 
  tm_raster (palette = "Purples", alpha = 0.6, style = "fixed", title = "% of area protected",
             breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
             labels = c(">0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100")) +
  tm_layout(legend.outside = TRUE, legend.outside.position = "right", outer.margins = c(0, 0, 0, 0))

GexpPD_hotspots.tm + ppa.tm + coastline.tm 

fig6 <- GexpPD_hotspots.tm + ppa.tm + coastline.tm 

png(filename = paste0(getwd(), "/outputs/fig6.png"),
    width = 13, height = 4, units = "cm", pointsize = 8, res = 600)
fig6
dev.off()

### Figure S1 - Correlations of species scores between the 2 phylogenies (PHYLACINE versus UPHAM) ----
scores_phylacine <- read.csv2(paste0(getwd(), "/outputs/scores_Mammals_median-over-the-100-resolved-trees.csv"))
scores_phylacine <- scores_phylacine[, c(1, 2, 4)]
colnames(scores_phylacine) <- c("species", "HEDGE_PHYLACINE", "LEDGE_PHYLACINE")
scores_upham <- read.csv2(paste0(getwd(), "/outputs/scores_Mammals_median-over-the-100-resolved-trees_UPHAM.csv"))
scores_upham <- scores_upham[, c(1, 2, 4)]
colnames(scores_upham) <- c("species", "HEDGE_UPHAM", "LEDGE_UPHAM")
scores_all <- merge(scores_phylacine, scores_upham, by = "species")

png(filename = paste0(getwd(), "/outputs/figS1a.png"),
    width = 17, height = 12, units = "cm", res = 600)
corPlot(df = scores_all[, c(2, 4)], method = "spearman", digits = 2, ties.method =  "average") # Spearman correlation for HEDGE
dev.off()

png(filename = paste0(getwd(), "/outputs/figS1b.png"),
    width = 17, height = 12, units = "cm", res = 600)
corPlot(df = scores_all[, c(3, 5)], method = "spearman", digits = 2, ties.method =  "average") # Spearman correlation for LEDGE
dev.off()

  
# comparing TOP HEDGE PHYLACINE and TOP HEDGE UPHAM
tophedge_upham <- head(scores_all[order(scores_all$HEDGE_UPHAM, decreasing = TRUE) ,], n = round(nrow(scores_all)*25/100))
species_tophedge_phylacine_upham <- intersect(tophedge$species_phylacine, tophedge_upham$species)
length(species_tophedge_phylacine_upham) # 1231 common species
length(species_tophedge_phylacine_upham)/nrow(tophedge_upham) *100 # 93 %
species_hedgeupham_nothedgephylacine <- setdiff(tophedge_upham$species, tophedge$species_phylacine)
species_hedgeupham_nothedgephylacine # 98 species TOP HEDGE UPHAM not in TOP HEDGE PHYLACINE
hedgeupham_nothedgephylacine <- tophedge_upham[tophedge_upham$species %in% species_hedgeupham_nothedgephylacine,]

# comparing TOP LEDGE PHYLACINE and TOP LEDGE UPHAM
topledge_upham <- head(scores_all[order(scores_all$LEDGE_UPHAM, decreasing = TRUE) ,], n = round(nrow(scores_all)*25/100))
species_topledge_phylacine_upham <- intersect(topledge$species_phylacine, topledge_upham$species)
length(species_topledge_phylacine_upham) # 736 common species
length(species_topledge_phylacine_upham)/nrow(topledge_upham) *100 # 55 %
species_ledgeupham_notledgephylacine <- setdiff(topledge_upham$species, topledge$species_phylacine)
species_ledgeupham_notledgephylacine # 593 species TOP HEDGE UPHAM not in TOP HEDGE PHYLACINE
ledgeupham_notledgephylacine <- topledge_upham[topledge_upham$species %in% species_ledgeupham_notledgephylacine,]
