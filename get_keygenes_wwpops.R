library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(stringr)


##############################################
#     DATA OCEANIA - FCS, files

path1 = "/gpfs/gibbs/pi/tucci/mpb63/infertility_2024/WWPops_bed_files"
path2 = "/gpfs/gibbs/pi/tucci/mpb63/OceanicPops_2024_FCS"

temp1 <- list.files(path = path1, pattern = "*0.0001.tsv", full.names = TRUE)
temp2 <- list.files(path = path2, pattern = "*0.0001.tsv", full.names = TRUE)

temp3 <- c(temp1, temp2)

data_list <- list()

#arch <- "Neandertal"

for (i in 1:length(temp3)) {
  name <- gsub(pattern = ".tsv", "", basename(temp3[i]))
  df <- read_tsv(temp3[[i]])
  df <- df %>%
    select(chr, start, end, population, FCS, log10_FCS, gene)
  df <- na.omit(df)
  #df <- df %>%
  # filter(str_detect(TractID, arch))
  data_list[[i]] <- df
}

merge_windows_within_10kb <- function(df) {
  df %>%
    arrange(chr, start) %>%
    group_by(chr) %>%
    mutate(
      overlapping = start < lag(end, default = first(end)),
      distance = start - lag(end, default = first(start)),
      group_id = cumsum(start > lag(end, default = first(start)) + 10000)
    ) %>%
    group_by(chr, group_id) %>%
    summarize(
      new_start = min(start),
      new_end = max(end),
      gene = paste(unique(unlist(strsplit(gene, ",", fixed =TRUE))),
                   sep=",", collapse = ","),
      population = first(population),
      FCS = first(FCS),
      log10_FCS = max(log10_FCS)
    ) %>%
    ungroup() %>%
    select(-group_id)
}

overlapping_10kb <- lapply(data_list, merge_windows_within_10kb)

# 10kb large windows per population
merged_overlapping_10kb <- bind_rows(overlapping_10kb)

###################################################################
# TOP WINDOWS PER POPULATION
###################################################################

# # GET the Top windows with the highest log10_FCS for each population
# top_windows_10kb_list <- lapply(overlapping_10kb, function(df) {
#   df %>%
#     arrange(desc(log10_FCS)) %>%
#     slice_head(n = 10) %>%
#     arrange(chr, new_start, new_end) %>%
#     ungroup()
# })
# 
# 
# # combine all populations top 20 windows per population
# merged_top_windows_10kb <- bind_rows(top_windows_10kb_list)

#view(merged_top_windows_10kb)

##############################################
##################################################

# KEY GENES SUBSET ## 
keygene_list <- read.table("infertilityGeneList_unique.txt")
keygene_list <- keygene_list[[1]]

match_genes <- function(gene_line) {
  genes <- strsplit(gene_line, ",")[[1]]
  any(genes %in% keygene_list)
}

# Subset keygene_list from overlapping_10kb
filtered_key_genes_10kb_list <- lapply(overlapping_10kb, function(df) {
  df %>%
    filter(sapply(gene, match_genes))
})

filtered_key_genes_10kb_list

###############
##############################################################################
##############################################################################
# IMPORTANT 
##################################################
# GET THE MERGED WINDOWS FOR THE OTHER POPULATIONS
##################################################

##################################################

#######################################################
# FUNCTION: 
path_ww = "/gpfs/gibbs/pi/tucci/mpb63/WWPops_2024_FCS"
path_oc = "/gpfs/gibbs/pi/tucci/mpb63/OceanicPops_2024_FCS"

tempoc <- list.files(path = path_oc, pattern = "*FCS.tsv", full.names = TRUE)
tempww <- list.files(path = path_ww, pattern = "*FCS.tsv", full.names = TRUE)

temp_all <- c(tempww, tempoc)

process_window <- function(full_df_path, interval_df) {
  
  full_df <- read_tsv(full_df_path, col_select = c(chr, start, end, population, FCS, log10_FCS, gene))
  #full_df <- read_tsv(temp[1], col_select = c(chr, start, end, population, FCS, log10_FCS, gene))
  # Create a new df
  #interval_df= bind_rows(filtered_key_genes_10kb_list)
  merged_df <- bind_rows(interval_df)
  # Set p0.0001 to TRUE in the merged_df data frame
  merged_df$p0.0001 <- TRUE
  
  # Create an empty data frame to store the merged rows and intervals
  merged_rows_df <- data.frame(chr = integer(),
                               new_start = integer(),
                               new_end = integer(),
                               population = character(),
                               gene = character(),
                               log10_FCS = numeric(),
                               p0.0001 = logical())
  
  # Filter and merge the rows from full_df that fall within the intervals specified in interval_df
  for (i in seq_len(nrow(interval_df))) {
    interval <- interval_df[i, ]
    merged_rows <- full_df %>%
      filter(chr == interval$chr, start >= interval$new_start, end <= interval$new_end) %>%
      mutate(gene = paste(unique(unlist(strsplit(gene, ",", fixed =TRUE))),
                          sep=",", collapse = ",")) %>%
      select(-c(start, end))  # Remove start and end columns
    
    # Add the interval values to merged rows
    merged_rows$new_start <- interval$new_start
    merged_rows$new_end <- interval$new_end
    
    # Append the current interval's merged rows to the result data frame
    merged_rows_df <- rbind(merged_rows_df, merged_rows)
  }
  
  # Combine the rows from interval_df and the merged rows with intervals
  merged_df <- bind_rows(merged_df, merged_rows_df)
  
  # Group by columns and summarize
  merged_df <- merged_df %>%
    group_by(chr, new_start, new_end, population, gene) %>%
    summarize(
      log10_FCS = max(log10_FCS, na.rm = TRUE),
      p0.0001 = all(p0.0001)
    ) %>%
    ungroup()
  
  merged_df$p0.0001[is.na(merged_df$p0.0001)] <- FALSE  
  
  rm(full_df);
  gc();
  
  return(merged_df)
}


#######################################################
# Apply function TO: 

#overlapping_10kb
#top_windows_10kb_list
#filtered_key_genes_10kb_list
###############################################
#######################################################

windows_key_genes_10kb <- lapply(temp_all, process_window, 
                                 interval_df= bind_rows(filtered_key_genes_10kb_list)) 

view(windows_key_genes_10kb[[1]])

merged_df <- bind_rows(windows_key_genes_10kb)

########################################################


####################################
#OPTION 1: ALL selective sweeps 

# merged_df <- merged_df[merged_df$p0.0001 == FALSE, ]
# 
# merged_df
# 

####################################
#OPTION 2: Arch selective sweeps 
merged_df <- bind_rows(windows_key_genes_10kb)

merged_df2 <- merged_df %>%
  group_by(chr, new_start, new_end, population, gene, log10_FCS) %>%
  arrange(desc(p0.0001)) %>%
  slice(1) %>%
  ungroup()

####################################
#OPTION 2: SECOND WINDOWED 

process_window2 <- function(df) {
  df %>%
    arrange(chr, new_start) %>%
    group_by(chr) %>%
    mutate(
      overlapping = new_start <= lag(new_end, default = first(new_start)),
      distance = new_start - lag(new_end, default = first(new_start)),
      group_id = cumsum(new_start >= lag(new_end, default = first(new_start)) & !overlapping)
    ) %>%
    group_by(chr,group_id, population) %>%
    summarize(
      new_start2 = min(new_start),
      new_end2 = max(new_end[length(new_end)]),
      log10_FCS = max(log10_FCS),
      gene = paste(unique(unlist(strsplit(gene, ",", fixed =TRUE))),
                   sep=",", collapse = ","),
      p0.0001_true = any(p0.0001),          # TRUE if at least one TRUE value in the group
    ) %>%
    ungroup() %>%
    mutate(
      new_coordinates2 = paste(chr, new_start2, sep = ":"),
      new_coordinates2 = paste(new_coordinates2, new_end2, sep = "-")
    ) %>%
    select(-group_id)
}

merged_df_key_genes <- process_window2(merged_df2)

#print(merged_df_key_genes)


############################################
# SAVE here and open in 03_positive_selection.Rmd
#############################################
save(merged_df, merged_df2, merged_df_key_genes, data_list, 
     file ="merged-windows2_10kb_heatmap_WWpop_FCS_0001.Rdata", 
     compress = "gzip", compression_level = 9)
