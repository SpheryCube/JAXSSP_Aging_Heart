# Author: Daniel Alfonsetti (daniel.alfonsetti@gmail.com)
# Date: July 10, 2018
# Description: This script reads output from the allez gsea output and QTL peak data.
# For each category in the gsea output, this script will look at the genes in that category
# and cross compare to see if those genes are in the corresponding QTL peak table.
# What I mean by "corresponding" is, for example, if the GSEA was done for protein expresison levels
# differentially expressed on age (denoted "protein_age") then the corresponding QTL peak table
# would be the pQTLs where 'age' was an interactive covariate (denoted "protein_age").

# Genes in the intersection of the two files are recorded and stored
# in another file. Using the genes in a given intersection, a QTL map is formed. 
# We are looking for sets of genes that share QTL peaks.
#######################
# Reset environment
rm(list = ls())
gc()
dev.off()

library(tidyverse) # Contains ggplot2, dplyr, stringr and some other stuff
library(devtools)
# Using local directories
data_file <- "~/do_heart/data/cleaned_data.RData"
load(data_file)
enrich_dir <- "/Users/c-alfond/do_heart/results/global_enrichment_allez/"
qtl_dir <- "/Users/c-alfond/do_heart/results/"
#######################
#######################
#######################
# Yuka Takemon and Dan Gatti's QTL scatterplot function
# Arguments:
# data: data.frame (or tibble) with the following columns:
#       ensembl: (required) character string containing the Ensembl gene ID.
#       qtl_chr: (required) character string containing QTL chromsome.
#       qtl_pos: (required) floating point number containing the QTL position 
#                in Mb.
#       qtl_lod: (optional) floating point number containing the LOD score.
#       gene_chr:  (optional) character string containing transcript chromosome.
#       gene_start: (optional) character string containing transcript start 
#                 postion in Mb.
#       gene_end:  (optional) character string containing transcript end
#                position in Mb.
# color.points: logical that is TRUE if the points should be colored by LOD.
# cis.points: logical that is TRUE if the points should be colored if they
#             are with in cis.
# cis.radius: numeric value containing the radius in Mb between a gene and a cis-eQTL.
#             Optional.
# cis.color: color for cis QTL. Optional.
# Returns:
# a plot of the QTL and gene location for each gene.

ggQTLmap = function(data, title, color.points = FALSE, cis.points = FALSE, cis.radius = 2, 
                    cis.color = "#4286f4") {
  
  
  # Check for required column names.
  required.colnames = c("ensembl", "qtl_chr", "qtl_pos")
  
  if(all(!required.colnames %in% colnames(data))) {
    stop(paste("colnames must contain the following columns:", 
               paste(required.colnames, collapse = ",")))
  } # if(!required.colnames %in% colnames(data))
  
  # Make sure that columns are not factors.
  data$ensembl = as.character(data$ensembl)
  data$qtl_chr = as.character(data$qtl_chr)
  
  gene.position.colnames = c("gene_chr", "gene_start", "gene_end")
  if(!all(gene.position.colnames %in% colnames(data))) {
    
    message(paste("Using Ensembl gene locations because optional gene",
                  "position columns (", paste0(gene.position.colnames, collapse = ","),
                  ") not found."))
    
    # Get the latest Ensembl GTF.
    ensembl = get_ensembl_genes()
    
    id    = ensembl$gene_id
    chr   = seqnames(ensembl)
    start = start(ensembl) * 1e-6
    end   = end(ensembl)   * 1e-6
    
    df = data.frame(ensembl = id, gene_chr = chr, gene_start = start, gene_end = end,
                    stringsAsFactors = F)
    data = left_join(data, df, by = "ensembl")
    
  } # if(gene.position.colnames %in% colnames(data))
  
  # Make sure that columns are not factors.
  data$gene_chr = as.character(data$gene_chr)
  
  # Get the gene mid-point.
  data = data %>% mutate(gene_pos = (gene_end + gene_start) * 0.5)
  
  # Fix the factor levels for the chr.
  all.chr = data %>% select(qtl_chr, gene_chr) %>%
    gather(k, v) %>%
    select(v) %>%
    distinct() %>%
    arrange(v)
  all.chr = all.chr$v[!is.na(all.chr$v)]
  
  if(length(grep("M", all.chr)) > 0) {
    wh = grep("M", all.chr)
    all.chr = all.chr[c(1:(wh-1), (wh+1):length(all.chr), wh)]
  }
  
  # Remove any NAs.
  data = data %>% na.omit
  
  data$qtl_chr  = factor(data$qtl_chr,  levels = all.chr[order(as.numeric(all.chr))])
  data$gene_chr = factor(data$gene_chr, levels = rev(all.chr[order(as.numeric(all.chr))]))
  
  # If we're plotting cis points, then add a cis-QTL column.
  if(cis.points) {
    
    data = data %>% mutate(cis = (gene_chr == qtl_chr) & (abs(gene_pos - qtl_pos) <= cis.radius))
    cis.colors = c("black", cis.color)
    names(cis.colors) = c("FALSE", "TRUE")
    out.plot <- ggplot(data, aes(x = qtl_pos, y = gene_pos)) +
      geom_point(aes(color = cis), alpha = 0.5) + 
      scale_color_manual(values = cis.colors) +
      facet_grid(gene_chr ~ qtl_chr, scales = "free", shrink = TRUE, drop = FALSE) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(fill = 0, color = "grey70"),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0.05, "lines"),
            axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(title)
    
  } else { 
    
    out.plot <- ggplot(data, aes(x = qtl_pos, y = gene_pos)) +
      geom_point(aes(color = qtl_lod, alpha = 0.5)) + {
        if(color.points) scale_color_continuous(low = "grey50", high = "red") 
      } +
      facet_grid(gene_chr ~ qtl_chr, scales = "free", shrink = TRUE, drop = FALSE) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(fill = 0, color = "grey70"),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0.05, "lines"),
            axis.text.x = element_text(angle = 90, hjust = 1)) + 
      ggtitle(title)
  } # else
  print(out.plot)
} # ggtmap()

#############################
#############################
# Part 1: Find genes that are in enriched categories and have QTLs for each test type.

setwd(qtl_dir)
peaks_files <- dir(pattern = "peaks")
setwd(enrich_dir)
gsea_files <- dir(pattern = "gsea_tbl_allez_clustered.csv")

test_types <- c("mrna_age", "mrna_sex", "protein_age", "protein_sex")

# Iterate over files
for (test_type in test_types) {
  
  peaks_file <- peaks_files[grepl(test_type, peaks_files)]
  gsea_file <- gsea_files[grepl(test_type, gsea_files)]
  
  # Read in data.
  setwd(qtl_dir)
  peaks_df <- read.csv(peaks_file)
  peaks_df <- peaks_df[peaks_df$qtl_lod > 6,]
  
  setwd(enrich_dir)
  gsea_tbl <- read.csv(gsea_file)
  
  
  summary_df <- data.frame(matrix(ncol = 3, nrow = nrow(gsea_tbl)))
  colnames(summary_df) <- c('group_id', 'genes_w_distal_qtl', 'genes_w_local_qtl')
  
  gsea_tbl_by_groups <- split(gsea_tbl, gsea_tbl$ClusterN)

  for (i in 1:length(gsea_tbl_by_groups)) {
    tbl <- gsea_tbl_by_groups[[i]]
    
    # Concatenate all the gene symbols from each category for this group
    ids <- paste(as.character(tbl$ensembl_ids), collapse = "")
    

    gsea_ensembl_ids <- strsplit(ids, ",")[[1]] 
    gsea_ensembl_ids <- gsea_ensembl_ids[!is.na(gsea_ensembl_ids)]
    
    
    peaks_df_row_indices <- if (substring(test_type, 1, 4) == "mrna") {
      match(gsea_ensembl_ids, peaks_df$lodcolumn) 
    } else {match(gsea_ensembl_ids, peaks_df$gene_id)} 
    peaks_df_row_indices <- peaks_df_row_indices[!is.na(peaks_df_row_indices)] # Remove NAs
    
    matched_rows <- peaks_df[peaks_df_row_indices,]
    # Split on distal and local qtl types.
    matched_rows_distal <- matched_rows[matched_rows$cis == FALSE,]
    matched_rows_local <- matched_rows[matched_rows$cis == TRUE,]
    matched_genes_distal <- as.character(matched_rows_distal$gene_symbol)
    matched_genes_local <- as.character(matched_rows_local$gene_symbol)
    
    # Put results into our summary df
    summary_df_row <- summary_df[i,]
    summary_df_row$group_id <- gsea_tbl_by_groups[[i]]$ClusterN[1]
    summary_df_row$genes_w_distal_qtl <- paste0(matched_genes_distal, collapse = ",")
    summary_df_row$genes_w_local_qtl <- paste0(matched_genes_local, collapse = ",")
    summary_df[i,] <- summary_df_row
  }
  setwd(enrich_dir)
  write.csv(summary_df, file = paste0("QTL_allezGSEA_cross_val_", test_type, ".csv"), row.names = FALSE)
}



#############################
#############################
# Part 2: Generate QTL maps againg, but now
# only plot the heatmap if a lot of genes have a peak on the same chromosome.
for (test_type in test_types) {
  print(paste0("Making QTLmaps for test type: ", test_type))
  # Set some parameters
  if (grepl("protein", test_type)) {expr_type = "protein"} else {expr_type = "mrna"}
  if (grepl("age", test_type)) {covar_name = "Age"} else if (grepl("sex", test_type)) {covar_name = "Sex"} else {covar_name = "None"}

  
  # Get peaks table
  setwd(qtl_dir)
  peaks_df <- read.csv(paste0("peaks_", test_type,"_thresh_6.csv")) 
  peaks_df <- peaks_df[peaks_df$qtl_lod > 6,]
  # Get genes that have a QTL and are in an enriched category (using data collected in part 1)
  setwd(enrich_dir)
  sig_qtls_by_group <- read.csv(paste0("QTL_allezGSEA_cross_val_", test_type, ".csv")) 
  
  setwd(enrich_dir)
  pdf(paste0("QTL_maps_for_", test_type,"_by_enriched_groups.pdf"))
  
  # Iterate over categories
  for (i in 1:nrow(sig_qtls_by_group)) {
    group_id <- sig_qtls_by_group[i,]$group_id # Category name
    
    # Get genes in current category.
    group_row <- sig_qtls_by_group[sig_qtls_by_group$group_id == group_id,]
    genes_1 <- strsplit(as.character(group_row$genes_w_distal_qtl), ",")[[1]]
    genes_2 <- strsplit(as.character(group_row$genes_w_local_qtl), ",")[[1]]
    genes <- c(genes_1, genes_2)
    # Pull out QTL peaks for genes in this category
    peaks_df_filtered <- peaks_df[peaks_df$gene_symbol %in% genes,] 
    
    if (nrow(peaks_df_filtered) == 0) {next}
    ###########################################################################
    # Create and save QTL map
    print("Saving QTL map...")
    if (expr_type == "mrna") { 
      title <- paste0( "eQTL Map (Interactive Covariate = ", covar_name, "):  Group ", i,  " Categories:")
      df <- rename(peaks_df_filtered, ensembl = lodcolumn)
    } else { 
      title <- paste0( "pQTL Map (Interactive Covariate = ", covar_name, "):  Group ", i, " Categories:")    
      df <- rename(peaks_df_filtered, ensembl = lodcolumn) 
    }
    ggQTLmap(data = df, title = title, color.points = FALSE, cis.points = TRUE, cis.radius = 2, 
             cis.color = "#4286f4")
    ###########################################################################
  } # for (i in 1:nrow(sig_qtls_by_group)) 
  dev.off()
  
} # for(test_type in test_types)