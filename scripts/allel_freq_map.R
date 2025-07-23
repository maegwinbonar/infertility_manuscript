## MAPS 

library(tidyverse)
library(sf)
theme_set(theme_bw())
library(ggspatial)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(mapview)
library(scatterpie)

######################################################
#Input: .csv file with the following columns: 
  #POP	- population
  #REF	- reference allele frequency
  #ALT  - alternate allele frequency
  #Latitude	-latitutde centroid from population samples
  #Longitude	- longitude centroid from population samples
  #rsid - rsid of lead variant
  #A1 - reference allele
  #A2 - alternate allele
  #effect - effect allele defined in GWAS
######################################################

global <- read.table("rsid_allele_freq.csv", sep=",", header = T)

world <- ne_countries(scale = "medium", returnclass = "sf", country = NULL)
class(world)

rsids <- unique(global$rsid)

for (i in rsids){
  df <- global %>%
    filter(rsid == i)
  
  ggplot(data = world) +
    geom_sf() +
    geom_scatterpie(data = df, aes(x = Longitude, y = Latitude, r = 7),
                    cols = c("REF", "ALT"), 
                    #alpha = 0.85,
                    colour = "black", linewidth=0.3,
                    col = NA) +
    scale_fill_manual(values = c("#FFFFA5", "#00008B"), 
                      labels = c(paste0("REF (", unique(df$A1), ")"), paste0("ALT (", unique(df$A2), ")"))) +
    ggrepel::geom_label_repel(data = df,
                             aes(x = Longitude, y = Latitude,
                                 label = POP),
                             colour = "black",
                             size = 4,
                             box.padding = 0.5,
                             segment.size = 0.5,
                             segment.color = 'black',
                             segment.alpha = 0.5,
                             force = 0.5,
                             min.segment.length = 0,
                             max.overlaps = Inf,
                             alpha = 0.8) +
    labs(fill = paste(df$rsid, df$effect, sep = ":")) +
    coord_sf(xlim = c(-150, 180), ylim = c(80, -60), expand = FALSE, 
             lims_method = c("cross"))+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "top")
  
  ggsave(paste0(i, "_rsid_allele_freq.pdf"), width = 8, height = 5)
}

ggplot(data = world) +
  geom_sf() +
  geom_scatterpie(data = global, aes(x = Longitude, y = Latitude),
                  cols = c("REF", "ALT"), 
                  alpha = 0.85, 
                  col = NA, 
                  stroke = 1) +
  scale_fill_manual(values = c("#99C24D", "#4BA5BE"), 
                    labels = c("Ancestral(A)", "Derived(G)")) +
  ggrepel::geom_text_repel(data = global,
                           aes(x = Longitude, y = Latitude,
                               label = POP),
                           colour = "black",
                           size = 3,
                           segment.size = 0.5,
                           segment.color = 'black',
                           segment.alpha = 0.5,
                           force = 0.5,
                           min.segment.length = 0,
                           max.overlaps = Inf) +
  labs(fill = "Allele") +
  coord_sf(xlim = c(-15, 180), ylim = c(80, -60), expand = TRUE, 
           lims_method = c("cross"))
