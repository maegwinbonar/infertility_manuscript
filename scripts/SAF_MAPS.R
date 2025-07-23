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
######################################################

oceania <- read.table("count_geo_Populations_combinedSAF_chr8-116372254.tsv", sep="\t", header = T)

oceania <- oceania %>%
  pivot_longer(
    cols = c(SAF_A, SAF_G),
    names_to = "SAF",
    values_to = "Percent"
  )

oceania

population_order <- c("Sepik","Goroka","Kove","Nakanai-Mangseng","Ata-Mamusi","Melamela","Baining-Kagat","Baining-Mali","Lavongai-Mussau",
                      "Nailik-Notsi-Tigak","Saposa","Nasioi","Vella-Lavella","Malaita","Bellona", "Rennell","Santa-Cruz","Tikopia")

pie_chart_list <- list()

for (group in population_order) {
  group_df <- oceania %>%
    filter(AnalysisGroup == group)
  
  pie_chart <- ggplot(group_df, aes(x = 1, y = Percent, fill = SAF)) +
    geom_bar(stat = "identity", width = 1, color = "grey20") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("SAF_A" = "#FEBF0F", "SAF_G" = "#F15B22")) +  
    labs(title = group) +
    theme_void() +
    theme(legend.position = "none")
  
  pie_chart_list[[group]] <- pie_chart
}

pie_chart_grid <- wrap_plots(pie_chart_list, nrow = 3)
pie_chart_grid



#######################################
# GLOBAL 
#######################################


global <- read.table("count_geo_Regions_chr8-116372254.tsv", sep="\t", header = T)

global <- global %>%
  pivot_longer(
    cols = c(SAF_A, SAF_G),
    names_to = "SAF",
    values_to = "Percent"
  )

print(global)

population_order <- c("AFR","AMR","CSA","EAS","EUR","ISEA","MDE","OCN")

pie_chart_list <- list()

for (group in population_order) {
  group_df <- global %>%
    filter(Region == group)
  
  pie_chart <- ggplot(group_df, aes(x = 1, y = Percent, fill = SAF)) +
    geom_bar(stat = "identity", width = 2, color = "grey20") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("SAF_A" = "#FEBF0F", "SAF_G" = "#F15B22")) +  
    labs(title = group) +
    theme_void() +
    theme(legend.position = "none")
  
  pie_chart_list[[group]] <- pie_chart
}

pie_chart_grid <- wrap_plots(pie_chart_list, nrow = 3)
pie_chart_grid



df_sf <- global %>%
  st_as_sf(coords = c("Latitude", "Longitude"), crs=4326)


world <- ne_countries(scale = "medium", returnclass = "sf", country = NULL)
class(world)


ggplot(data = world) +
  geom_sf() +
  geom_point(data = global, aes(x = Longitude, y = Latitude, 
                                fill=Region),
             position=position_jitter(h=0.1, w=0.1),
             shape = 21, alpha = 0.7, size = 3) +
  ggrepel::geom_label_repel(data = global, 
                            aes(x = Longitude, y = Latitude,
                                label = Region), 
                            colour = "black",
                            size = 2,
                            #segment.colour = "darkgrey",
                            min.segment.length = 0, 
                            max.overlaps = 20) +
  coord_sf(xlim = c(-125, 170), ylim = c(-50, 88), expand = TRUE, 
           lims_method = c("cross"))


ggplot(data = world) +
  geom_sf() +
  geom_scatterpie(data = global, aes(x = Longitude, y = Latitude),
                  cols = c("SAF_A", "SAF_G"), 
                  alpha = 0.85, 
                  col = NA, 
                  stroke = 1) +
  scale_fill_manual(values = c("#99C24D", "#4BA5BE"), 
                    labels = c("Ancestral(A)", "Derived(G)")) +
  ggrepel::geom_text_repel(data = global,
                           aes(x = Longitude, y = Latitude,
                               label = Region),
                           colour = "black",
                           size = 3,
                           segment.size = 0.5,
                           segment.color = 'black',
                           segment.alpha = 0.5,
                           force = 0.5,
                           min.segment.length = 0,
                           max.overlaps = Inf) +
  labs(fill = "Allele") +
  coord_sf(xlim = c(-125, 170), ylim = c(-50, 88), expand = TRUE, 
           lims_method = c("cross"))

######################################################
# Maegwin's Data
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
