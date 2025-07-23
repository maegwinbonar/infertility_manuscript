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
library("EnsDb.Hsapiens.v75")
db <- EnsDb.Hsapiens.v75

#ensDb_GRCh37_73 <- ah[["AH10684"]]


path = "/gpfs/gibbs/pi/tucci/mpb63/WWPops_2024_FCS"
pattern1 <- "FCS.tsv"
#genes <- read.table("infertility_genes.txt", header = TRUE)

# list of files
temp <- list.files(path = path, pattern = pattern1, full.names = TRUE)

#data_list <- list()

# gene list
genes <- read_tsv("infertility_genes_FCS_0.0001.tsv", col_names = FALSE)
#genes$pos <- as.numeric(genes$pos)
gene_list <- unique(genes[[1]])
# metadata
genes_df <- read_delim("infertilityGeneList.txt")


for (i in 1:length(temp)) {
#name <- gsub(pattern = ".tsv", "_locuszoom", basename(temp[i]))
  name_plot <- gsub(pattern = "_PIBv1.*", "", basename(temp[i]))

  df <- read_tsv(temp[i]) %>%
  dplyr::select(chr, start, end, avg_pos, population, FCS, log10_FCS, gene) %>%
  na.omit()

  #data_list[[i]] <- df

# make snp id
  df$SNP <- paste(df$chr, df$avg_pos, sep = "_")

# make avg position a whole number
  df$pos <- floor(df$avg_pos)

#Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('dodgerblue','green'))

#This adds a column of color values
# based on the y values
  df$Col <- rbPal(4)[as.numeric(cut(df$log10_FCS,breaks = 4))]

for (j in 1:length(gene_list)) {

# locus zoom function
  loc <- locus(data = as.data.frame(df), 
             gene = gene_list[j], 
             flank = 10e5,
             ens_db = db,
             chrom = 'chr',
             pos = 'pos',
             labs = 'SNP',
             yvar = 'log10_FCS')


# plot it

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
    dplyr::filter(chr %in% dt$chr) %>%
    dplyr::filter(`pos(GRCh37)` > min(dt$start) & `pos(GRCh37)` < max(dt$end))
  
  # binary search and "roll" to the nearest neighbour
  # In the final expression the val column will have the you're looking for.
  var <- dt[J(as.numeric(var_highlight$`pos(GRCh37)`)), roll = "nearest"]
  
  var <- var %>%
    left_join(var_highlight, by=c("val"="pos(GRCh37)"))
  
  if (nrow(var) > 0){
    pl <- quote({
      points(x=var$pos, y=var$log10_FCS, pch = 25, bg = "red")
      # add custom text label for index SNP
      lx <- var$pos
      ly <- jitter(var$log10_FCS+3, amount = 0.5)
      text(lx, ly, var$rsid, pos = 4, cex = 0.8)
      arrows(var$pos, var$log10_FCS, lx+3e4, ly, length = 0, lty = "dashed")
    })
    
    pdf(file = paste0("plots/locus_zoom/", gene_list[j], "/", name_plot, "_", gene_list[j], ".pdf"), width = 12, height = 8, pointsize = 15)
    locus_plot(loc, pch = 21,
               bg=loc$data$Col,
               index_snp = NULL,
               pcutoff = NULL, panel.first = pf, panel.last=pl,
               highlight = var$gene.x,
               highlight_col = "red3",
               ylim = c(0,6))
    title(main = paste(name_plot, "-", gene_list[j], sep = " "), adj = 0, line = 2.5)
    dev.off()
    
  } else {
    # plot using locus_plot
    # save it
    pdf(file = paste0("plots/locus_zoom/", gene_list[j], "/", name_plot, "_", gene_list[j], ".pdf"), width = 12, height = 8, pointsize = 15)
    locus_plot(loc, pch = 21,
               bg=loc$data$Col,
               index_snp = NULL,
               pcutoff = NULL, panel.first = pf,
               highlight = gene_list[[j]],
               highlight_col = "red3",
               ylim = c(0,6))
    title(main = paste(name_plot, "-", gene_list[j], sep = " "), adj = 0, line = 2.5)
    
    dev.off()
  }
  
}

}

# 