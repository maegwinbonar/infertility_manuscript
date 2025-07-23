##############################################
#     DATA World Wide Populations - locuszoom plots of the 21 infertility genes plot | FCS

#devtools::install_github("myles-lewis/locuszoomr")

library(locuszoomr)
library(dplyr)
library(readr)
library(data.table)
library(tidyverse)
library(stringr)
library(AnnotationHub)
library(ensembldb)

#ah <- AnnotationHub()
#query(ah, c("EnsDb.Hsapiens.v75"))
#ensDb_GRCh37_73 <- ah[["AH10684"]]

## Annotation from ensembldb to get one with compatible with hg37
library("EnsDb.Hsapiens.v75")
db <- EnsDb.Hsapiens.v75

## Path to FCS files
path1 = "/gpfs/gibbs/pi/tucci/mpb63/OceanicPops_2024_FCS"
path2 = "/gpfs/gibbs/pi/tucci/mpb63/WWPops_2024_FCS"

# gene list
genes <- read_table("r2_selection_windows.txt", col_names = FALSE)
colnames(genes) <- c("gene", "pop", "chr", "start", 'end', 'var')
#genes$pos <- as.numeric(genes$pos)
gene_list <- unique(genes[[1]])
# metadata
genes_df <- read_delim("infertilityGeneList.txt")

## Looped over each population
for (i in 1:nrow(genes)) {
  
  # make a path file with the correct path
  if(grepl("OCN", genes$pop[i])){
    temp <- list.files(path = path1, pattern = paste0(genes$pop[i], "_", ".*FCS.tsv$$"), full.names = TRUE)
 } else{
    temp <- list.files(path = path2, pattern = paste0(genes$pop[i], "_", ".*FCS.tsv$"), full.names = TRUE)
 }
  
  #name <- gsub(pattern = ".tsv", "_locuszoom", basename(temp[i]))
  name_plot <- gsub(pattern = "_PIBv1.*", "", basename(temp))
  
  df <- read_tsv(temp) %>%
    dplyr::select(chr, start, end, avg_pos, FCS, log10_FCS, gene) %>%
    dplyr::filter(chr == genes$chr[i]) %>%
    na.omit()
  
  # make snp id
  df$SNP <- paste(df$chr, df$avg_pos, sep = "_")
  
  # make avg position a whole number
  df$pos <- floor(df$avg_pos)
  
  #Get LD values
  ld_tmp <- list.files(path = "r2", pattern = paste0(genes$pop[i], "_", genes$gene[i], ".ld"), full.names = TRUE) 
  ld_df <- read.table(ld_tmp, header = TRUE)
  
  setDT(ld_df)
  setDT(df)
  #make r2 clumn
  df[, r2 := NA_real_]
  
  # add r2 values
  for (j in seq_len(nrow(ld_df))) {
    bp <- ld_df[j, BP_B]
    r2_val <- ld_df[j, R2]
    
    # Find matching rows
    idx <- which(df$start <= bp & df$end >= bp)
    
    # For each match, update with max r2 value
    for (k in idx) {
      current_r2 <- df[k, r2]
      if (is.na(current_r2) || r2_val > current_r2) {
        df[k, r2 := r2_val]
      }
    }
  }
  
  ## for each population loop over each gene
  # locus zoom function
    
  if(genes$gene[i]=="ZEB2"){
    loc <- locus(data = as.data.frame(df),
                 gene = genes$gene[i], 
                 flank = 10e5,
                 ens_db = db,
                 chrom = 'chr',
                 pos = 'pos',
                 labs = 'SNP',
                 yvar = 'log10_FCS',
                 LD='r2')
  }else{
  loc <- locus(data = as.data.frame(df),
                 gene = genes$gene[i], 
                 flank = 5e5,
                 ens_db = db,
                 chrom = 'chr',
                 pos = 'pos',
                 labs = 'SNP',
                 yvar = 'log10_FCS',
                 LD='r2')
    
  } 
    
    # create a panel with the dashed line at significance value 2 to add to plot
    pf <- quote({
      abline(h = 4, col = "grey50", lty = 2)
      abline(h = 5, col = "grey50", lty = 2)
    })
    
    ### this is to highlight my specific rsids of interest
    ### finds the windows which encompasses the rsid
    dt = data.table(loc$data, val = loc$data$pos) # you'll see why val is needed in a sec
    setattr(dt, "sorted", "x")  # let data.table know that w is sorted
    setkey(dt, val) # sorts the data
    
    # create a point with the variant from the GWAS highlighted
    
    var_highlight <- genes_df %>%
      dplyr::filter(rsid == genes$var[i]) %>%
      dplyr::filter(`pos(GRCh37)` > min(dt$start) & `pos(GRCh37)` < max(dt$end))
    
    # binary search and "roll" to the nearest neighbour
    # In the final expression the val column will have the you're looking for.
    var <- dt[J(as.numeric(var_highlight$`pos(GRCh37)`)), roll = "nearest"]
    
    var <- var %>%
      left_join(var_highlight, by=c("val"="pos(GRCh37)"))
    
      pl <- quote({
        points(x=var$pos, y=var$log10_FCS, pch = 23, bg = "yellow", cex = 1.1)
        # add custom text label for index SNP
        lx <- var$pos
        ly <- jitter(var$log10_FCS+3, amount = 0.5)
        text(lx, ly, var$rsid, pos = 4, cex = 0.8)
        arrows(var$pos, var$log10_FCS, lx+3e4, ly, length = 0, lty = "dashed")
      })
      
      ## save as png
      pdf(file = paste0("plots/locus_zoom/r2/", name_plot, "_", genes$gene[i], "_r2.pdf"), width = 12, height = 8, pointsize = 15)
      locus_plot(loc, pch = 21,
                 index_snp = NULL,
                 # See add the first and last panels
                 pcutoff = NULL, panel.first = pf, panel.last=pl,
                 highlight = var$gene.x,
                 highlight_col = "red3",
                 ylim = c(0,6))
      title(main = paste(name_plot, "-", genes$gene[i], sep = " "), adj = 0, line = 2.5)
      
      dev.off()

    
    
    
}


