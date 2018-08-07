# Find and analyze eQTL and pQTL hotspots for DO cross sectional heart data.
# Create QTL density plots and QTL maps then perform GSEA with the genes that map to the hotspots.
# R/3.4.4
# Daniel Alfonsetti
# daniel.alfonsetti@gmail.com
# July. 30, 2018
################################################################################
#######################
# Reset environment
rm(list = ls())
gc()

library(qtl2)
library(broom)
library(RColorBrewer)
library(doParallel)
library(pcaMethods)

library("org.Mm.eg.db")  # Mouse data base
library(tidyverse) # Contains ggplot2, dplyr, stringr and some other stuff

# Enrichment analysis libraries
# library(allez)
library(clusterProfiler)

######## Local directories
num_cores <- 2

# Input paths
peaks_tables_dir <- "~/do_heart/results/"
perms_input_dir <- "~/do_heart/results/"

# Output paths
base_output_dir <- "~/do_heart/results/hotspots/"
dir.create(base_output_dir)
load("~/do_heart/data/cleaned_data.RData")

######## HPC directories
# num_cores <- 30
# #
# # Input paths
# peaks_tables_dir <- "/fastscratch/c-alfond/do_heart/results/"
# perms_input_dir <- "/fastscratch/c-alfond/do_heart/results/"
# #
# # # Output paths
# base_output_dir <- "/fastscratch/c-alfond/do_heart/results/hotspots/"
# dir.create(base_output_dir)
# load("/fastscratch/c-alfond/do_heart/data/cleaned_data.RData")

##########################################################
# Register parallele back end.
registerDoParallel(cores = num_cores)
##########################################################
# Helper functions

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

# Dan Gatti's QTL Density Plot function
qtl_density_func <- function(peaks_table, lod_threshold, cis_boolean, expr_type, covar_name, test_type, create_file = TRUE) {
  
  # Window 4Mb, slide it in 1Mb steps.
  breaks = matrix(c(seq(0, 200, 4), seq(1, 201, 4), seq(2, 202, 4), seq(3, 203, 4)), ncol = 4)
  tmp = as.list(1:ncol(breaks)) 
  for(i in 1:ncol(breaks)) {
    tmp[[i]] = peaks_table %>%
      filter(qtl_lod >= lod_threshold & cis == cis_boolean) %>% 
      arrange(qtl_chr, qtl_pos) %>%
      group_by(qtl_chr) %>%
      mutate(win = cut(qtl_pos, breaks = breaks[,i])) %>%
      group_by(qtl_chr, win) %>% 
      summarize(cnt = n()) %>%
      separate(win, into =  c("placeholder", "prox", "dist")) %>%
      mutate(prox = as.numeric(prox), 
             dist = as.numeric(dist), 
             mid = 0.5 * (prox + dist)) %>%
      select(qtl_chr, mid, cnt)
  }
  
  result = bind_rows(tmp[[1]], tmp[[1]], tmp[[3]], tmp[[4]])
  rm(tmp)
  
  if (cis_boolean & expr_type == "mrna") {
    plot_title <- paste0("cis-eQTL Density Histogram (Interactive Covariate=", covar_name, ")"); 
    file_name <- "cis_eqtl_density_hist.pdf"}
  else if (cis_boolean & expr_type == "protein") {
    plot_title <- paste0("cis-pQTL Density Histogram (Interactive Covariate=", covar_name, ")"); 
    file_name <- "cis_pqtl_density_hist.pdf"}
  else if(!cis_boolean & expr_type == "mrna") {
    plot_title <- paste0("trans-eQTL Density Histogram (Interactive Covariate=", covar_name, ")"); 
    file_name <- "trans_eqtl_density_hist.pdf"}
  else {
    plot_title <- paste0("trans-pQTL Density Histogram (Interactive Covariate=", covar_name, ")"); 
    file_name <- "trans_pqtl_density_hist.pdf"}
  
  # Make plot
  out.plot = ggplot(result, aes(mid, cnt, color = cnt)) +
    geom_line() +
    scale_color_gradient(low = "#32b72d", high = "#af33c6") + 
    facet_grid(.~qtl_chr, scales = "free") +
    theme_light() +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = 0, color = "grey70"),
          panel.spacing = unit(0, "lines"),
          axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title = plot_title, x = "Mb", y = "Number of QTLs")
  
  if (create_file) {
    file_name <- paste0(test_type,"_", file_name)
    pdf(paste0(output_dir, file_name), width = 10, height = 8)
    print(out.plot)
    dev.off()
  } else {
    print(out.plot)
  }
  return(result)
}


# qtl_type: is in reference to the qtl location. 'cis' vs 'trans'
# scan_type: is in reference to the model and expression data used. 
# For example, 'protein_age' means the QTLs are pQTLs where age was an interactive covariate.
# expr_to_plot: Either 'mrna' or 'protein'. If we have eQTL data and we want to plot protein,
# the function will plot proteins that correspond to the genes in the eQTL data. Vice versa
# if we have pQTL data and we want to plot their transcript correlations.

hotspot_correlations <- function(hotspot.genes, qtl_type, expr_to_plot = "mrna", test_type){

  # Hotspot correlations
  breaks = -100:100/100
  colors = colorRampPalette(rev(brewer.pal(11, "Spectral")))(length(breaks) - 1)

  file_name <-  paste0(test_type, "_", qtl_type, "_qtl_gene_hotspot_", expr_to_plot, "_correlations.pdf")
  pdf(paste0(output_dir, file_name), width = 10, height = 10)

  for(i in 1:length(hotspot.genes)) {
    chr = names(hotspot.genes)[i]

    if (length(hotspot.genes[[i]]$lodcolumn) == 1) {next} # Can't do correlations if there is only one gene in hotspot.

    if (expr_to_plot == "protein") { # plot protein-protein correlations for genes in the hotspot

      # If the lod peaks are for genes and we want to plot the correlations between those genes' corresponding proteins
      if (grepl("mrna", test_type)) {

        # Get the corresponding protein ids and then plot. For each protein id, record the qtl score
        # of the gene that mapped to it so that we can color the protein ids on the heatmap.
        # Note: we don't have protein data that corresponds to every mrna. Therefore we will just skip
        # those that don't have corresponding mrna.
        df <- data.frame(matrix(nrow = 0, ncol = 2))
        colnames(df) <- c("protein_id", "corresponding_gene_qtl")
        for (j in 1:nrow(hotspot.genes[[i]])) { # for each gene...
          cur_gene_id <- hotspot.genes[[i]][j,]$lodcolumn
          cur_gene_qtlscore <- hotspot.genes[[i]][j,]$qtl_lod
          protein_ids <- annot.protein %>% filter(gene_id == cur_gene_id)  %>% select(protein_id)
          if (nrow(protein_ids) == 0) {next}
          for (k in 1:nrow(protein_ids)) # For each protein that this gene maps to...
          {
            protein_id <- protein_ids$protein_id[k]
            df[nrow(df)+1,] <- c(protein_id, cur_gene_qtlscore)
          } #for(k)
        } #for(j)
        df$corresponding_gene_qtl <- as.numeric(df$corresponding_gene_qtl)
        df <- df[!duplicated(df$protein_id),] # Get rid of repeats. There might be repeats since some genes can map to the same protein.
        if (nrow(df) == 0) {next} # Make sure we aren't empty now.
        side.colors = cut(df$corresponding_gene_qtl, breaks = 100)

        # make correlations matrix now.

        # Not all gene ids have corresponding proteins, so only get the genes that do.
        protein_ids <- annot.protein %>% filter(gene_id %in% hotspot.genes[[i]]$lodcolumn) %>% select(protein_id)
        protein_ids <- unlist(protein_ids)
        expr.data <- expr.protein[, protein_ids]
        if(length(expr.data) == 189) {next} # If we only have one viable protein, we can't plot anything.
        tmp <- cor(expr.data,  use = "pairwise.complete.obs")
        dimnames(tmp) = list(protein_ids, protein_ids)

        # If the lod peaks are for proteins and we want to plot protein correlation...
      } else {

        expr.data <- expr.protein[, hotspot.genes[[i]]$lodcolumn]
        tmp <- cor(expr.data,  use = "pairwise.complete.obs")
        dimnames(tmp) = list(hotspot.genes[[i]]$lodcolumn, hotspot.genes[[i]]$lodcolumn)
        side.colors = cut(hotspot.genes[[i]]$qtl_lod, breaks = 100)

      }

    } else { # We want to plot mrna - mrna correlations for genes in the hotspot

      # If the lod peaks are for genes and we want to plot gene correlations...
      if (grepl("mrna", test_type)) {

        expr.data <- expr.mrna[, hotspot.genes[[i]]$lodcolumn]
        tmp <- cor(expr.data,  use = "pairwise.complete.obs")
        dimnames(tmp) = list(hotspot.genes[[i]]$lodcolumn, hotspot.genes[[i]]$lodcolumn)
        # side.colors = cut(hotspot.genes[[i]]$qtl_lod, breaks = 100)
        side.colors = cut(hotspot.genes[[i]]$qtl_lod, breaks = 100)

        # If the lod peaks are for proteins and we want to plot mrna-mrna correlations...
      } else {

        # Get the corresponding gene ids and then plot
        expr.data <- expr.mrna[, hotspot.genes[[i]]$gene_id]
        tmp <- cor(expr.data,  use = "pairwise.complete.obs")
        dimnames(tmp) = list(hotspot.genes[[i]]$gene_id, hotspot.genes[[i]]$gene_id)
        side.colors = cut(hotspot.genes[[i]]$qtl_lod, breaks = 100)
      }
    }
    
    ########
    # Make heatmap
    ########

    # Make heatmap title
    if (grepl("age", test_type)) {covar_name = "Age"} else if (grepl("sex", test_type)) {covar_name = "Sex"} else {covar_name = "None"}
    if (grepl("mrna", test_type)) {
      plot_title <- paste0(qtl_type, " eQTLs Chr", chr, " Gene Hotspot ", expr_to_plot, " Correlations (Interactive Covariate = ", covar_name, ")")
    } else {
      plot_title <- paste0(qtl_type, " pQTLs Chr", chr, " Gene Hotspot ", expr_to_plot, " Correlations (Interactive Covariate = ", covar_name, ")")
    }

    side.colors = colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(length(levels(side.colors)))[as.numeric(side.colors)]
    names(side.colors) = rownames(tmp)

    # Plot
    heatmap(tmp, symm = TRUE, scale = "none", main = plot_title,
            breaks = breaks, col = colors, RowSideColors = side.colors, ColSideColors = side.colors)
  }
  dev.off()
}


get_hotspots_genes_func <- function(peaks_table, hotspots, qtl_type, lod_threshold, test_type) {
  
  hotspots_genes = as.list(hotspots$qtl_chr)
  names(hotspots_genes) = hotspots$qtl_chr
  
  if (grepl("protein", test_type)) {expr.data <- expr.protein
  } else {expr.data <- expr.mrna}
  
  
  for(i in 1:nrow(hotspots)) {
    
    if (qtl_type == "trans") {
      hotspots_genes[[i]] = peaks_table %>% 
        filter(qtl_lod >= lod_threshold) %>%
        filter(qtl_chr == hotspots$qtl_chr[i] & 
                 qtl_pos >= hotspots$proximal[i] & 
                 qtl_pos <= hotspots$distal[i]) %>%
        filter(cis == FALSE)
    } else {
      hotspots_genes[[i]] = peaks_table %>% 
        filter(qtl_lod >= lod_threshold) %>%
        filter(qtl_chr == hotspots$qtl_chr[i] & 
                 qtl_pos >= hotspots$proximal[i] & 
                 qtl_pos <= hotspots$distal[i]) %>%
        filter(cis == TRUE)
    }
    
    
    if (nrow(hotspots_genes[[i]]) == 0) {next}
    # Append current gene's expression data to the row.
    hotspots_genes[[i]][, rownames(expr.data)] <- NA
    for (j in 1:nrow(hotspots_genes[[i]])) {
      ensembl_id <- hotspots_genes[[i]]$lodcolumn[j]
      hotspots_genes[[i]][j, rownames(expr.data)] <- expr.data[,as.character(ensembl_id)]
    } # for(i)
    write_csv(hotspots_genes[[i]], path = paste0(output_dir, test_type, "_", qtl_type, "_chr", names(hotspots_genes)[i], "_hotspot_genes.csv"))
  } # for(j)
  return(hotspots_genes)
} #get_hotspots_genes_func


allez_helper_func <- function(allez_table) {
  # Description: Takes a table from the output of the allez function and adds columns that contain gene symbols and ensembl ids 
  # corresponding to the set of entrenz gene ids in a given row (category)
  ####
  
  # Convert entrenz gene ids into ensemble gene ids and gene symbols
  allez_table[, c("ensembl_ids", "symbols")] <- NA
  names(allez_table)[names(allez_table) == "genes"] <- "entrenz_ids"
  
  # Iterate over each row (category) and update them (add gene ids and gene symbols)
  for (i in 1:nrow(allez_table)) {
    # i = 1
    row <- allez_table[i,]
    row$entrenz_ids <- gsub(";", ",", row$entrenz_ids)
    
    entrenz_ids <- strsplit(row$entrenz_ids, ",")[[1]] 
    
    #result <- select(org.Mm.eg.db, entrenz_ids, c("SYMBOL", "GENENAME"), "ALIAS")
    result <- AnnotationDbi::select(org.Mm.eg.db, 
                                    keys= as.character(as.matrix(entrenz_ids)), 
                                    columns=c("SYMBOL", "ENSEMBL"), 
                                    keytype="ENTREZID")
    result <- result[!is.na(result$ENTREZID),]
    
    row$ensembl_ids <-  paste0(result$ENSEMBL, collapse = ",")
    row$symbols <- paste0(result$SYMBOL, collapse = ",")
    allez_table[i,] <- row
  }
  
  allez_table <- allez_table[order(allez_table[,'z.score'], decreasing = TRUE),] # Sort by z-score
  return(allez_table)
}

# End of helper functions
#####################
#####################
#########################################################################################################
#########################################################################################################
#########################################################################################################


output_dir <- paste0(base_output_dir, "preliminary_hotspot_results/")
dir.create(output_dir)
setwd(output_dir)

# Iterate over each test type
peaks_files <- Sys.glob(paste0(peaks_tables_dir, "*peaks*thresh_6*")) # contains the full path
peaks_files <- peaks_files[!grepl("full", peaks_files)]
print(peaks_files)

pdf("preliminary_hotspot_results.pdf")
# foreach (test_type = c("mrna_none", "mrna_sex", "mrna_age", "protein_none", "protein_sex", "protein_age")) %dopar% {
for (test_type in c("mrna_none", "mrna_sex", "mrna_age", "protein_none", "protein_sex", "protein_age"))  {
  
  # Set some parameters
  if (grepl("protein", test_type)) {expr_type = "protein"} else {expr_type = "mrna"}
  if (grepl("age", test_type)) {covar_name = "Age"} else if (grepl("sex", test_type)) {covar_name = "Sex"} else {covar_name = "None"}
  
  
  # Read in peaks table
  print(test_type)
  file <- peaks_files[grepl(test_type, peaks_files)]
  print(file)
  peaks_table <- read.csv(file, stringsAsFactors = TRUE); print("Read in table!")
  
  levels = c(1:19, "X", "Y", "M")
  peaks_table$gene_chr <- factor(peaks_table$gene_chr, levels = levels)
  peaks_table$qtl_chr <- factor(peaks_table$qtl_chr, levels = levels)
  
  # If using permutation testing...
    # Filter out genes with qtl's less than the 95 % of distribution maximum LODs
    # of 10000 permutations
    # Get current test types' corresponding set of permutations.
    # perm_file_name <- paste0("perms_", test_type, ".rds")
    # perms <- readRDS(paste0(perms_input_dir, perm_file_name))
    # 
    # lod_threshold <- summary(perms)[,'pheno1'] # 95 percentile
  # If using flat cutoff...
    lod_threshold <- 7.2
  peaks_table <- peaks_table %>% filter(qtl_lod > lod_threshold)
  
  ###########################################################################
  # Create and save QTL map
  print("Saving QTL map...")
  if (expr_type == "mrna") { 
    df <- rename(peaks_table, ensembl = lodcolumn)
    title <- paste0("eQTL Map (Interactive Covariate = ", covar_name, ")")
  } else { 
    title <- paste0("pQTL Map (Interactive Covariate = ", covar_name, ")")
    df <- rename(peaks_table,ensembl = lodcolumn) 
  }
  ggQTLmap(data = df, title = title, color.points = FALSE, cis.points = TRUE, cis.radius = 2, 
           cis.color = "#4286f4")
  ###########################################################################
  # Create and save density histograms
  print("Creating density histograms...")
  for (qtl_type in c("cis", "trans")) {
    cis_boolean <- grepl("cis", qtl_type)
    qtl_density_result <- qtl_density_func(peaks_table, lod_threshold = lod_threshold, cis_boolean = cis_boolean, 
                                           expr_type = expr_type, covar_name = covar_name, test_type = test_type, create_file = FALSE)
  }
}
dev.off()

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
# Manually set threshold for number of counts in a hotspot for each test type based on qtl density histograms
# made above. Each test vector is of the following form:
# c("TEST_TYPE_NAME", "CIS_HOTSPOT_CNT_THRESHOLD", "TRANS_HOTSPOT_CNT_THRESHOLD")


tests <- list(c("mrna_none", 80, 75), c("mrna_sex", 0, 40), c("mrna_age", 0, 25),
  c("protein_none", 0, 12), c("protein_sex", 0, 15), c("protein_age", 0, 35))

# Iterate over each test. For each hotspot get hotspot genes,
# make principal component plots, and do enrichment analysis

for (test_params in tests) {
# for (test_params in tests) {
  
  print(test_params)
  test_type <- test_params[1]
  if (grepl("protein", test_type)) {expr_type = "protein"} else {expr_type = "mrna"}
  if (grepl("age", test_type)) {covar_name = "Age"} else if (grepl("sex", test_type)) {covar_name = "Sex"} else {covar_name = "None"}
  
  
  file <- peaks_files[grepl(test_type, peaks_files)]
  print(file)
  dir.create(paste0(base_output_dir, test_type, "_hotspots_output/"))
  output_dir <- paste0(base_output_dir, test_type, "_hotspots_output/")
  setwd(output_dir)
  getwd()
  
  # Read in peaks table
  print(file)
  peaks_table <- read.csv(file, stringsAsFactors = TRUE)
  levels = c(1:19, "X", "Y", "M")
  print("Got here")
  peaks_table$gene_chr <- factor(peaks_table$gene_chr, levels = levels)
  print("Didn't get here")
  peaks_table$qtl_chr <- factor(peaks_table$qtl_chr, levels = levels)
  
  # If using permutations...
    # Filter out genes with qtl's less than the 95 % of distribution maximum LODs
    # of 10000 permutations
    # Get current test types' corresponding set of permutations.
    # perm_file_name <- paste0("perms_", test_type, ".rds")
    # perms <- readRDS(paste0(perms_input_dir, perm_file_name))
    # lod_threshold <- summary(perms)[,'pheno1'] # 95 percentile
  # If using flat cutoff...
    lod_threshold <- 7.2
  
  peaks_table <- peaks_table %>% filter(qtl_lod > lod_threshold)
  
  ###########################################################################
  # Create and save QTL map
  print("Making QTL map...")
  pdf(paste0(test_type, "_QTLmap.pdf"))
  if (expr_type == "mrna") { 
    df <- rename(peaks_table, ensembl = lodcolumn)
    title <- paste0("eQTL Map (Interactive Covariate = ", covar_name, ")")
  } else { 
    title <- paste0("pQTL Map (Interactive Covariate = ", covar_name, ")")
    df <- rename(peaks_table,ensembl = lodcolumn) 
  }
  ggQTLmap(data = df, title = title, color.points = FALSE, cis.points = TRUE, cis.radius = 2, 
           cis.color = "#4286f4")
  dev.off()
  ###########################################################################
  
  # Stratify by QTL position type.
  for (qtl_type in c("trans")) {
    print(paste0("QTL type ", qtl_type))
    
    
    ###########################################################################
    # Create and save density histograms
    cis_boolean <- grepl("cis", qtl_type)
    qtl_density_result <- qtl_density_func(peaks_table, lod_threshold = lod_threshold, cis_boolean = cis_boolean, 
                                           expr_type = expr_type, covar_name = covar_name, test_type = test_type, create_file = TRUE)
    ###########################################################################
    print("Saved density histograms...")
    
    # Pull out hotspot count threshold for this type.
    if (grepl("cis", qtl_type)) {hotspot_cnt_thresh <- as.integer(test_params[2])} else
    {hotspot_cnt_thresh <- as.integer(test_params[3])}
    
    if (grepl("protein", test_type)) {expr_type = "protein"} else {expr_type = "mrna"}
    if (grepl("age", test_type)) {covar_name = "Age"} else if (grepl("sex", test_type)) {covar_name = "Sex"} else {covar_name = "None"}
    
    ###########################################################################
    ###########################################################################
    # Get hotspots (i.e. genome intervals of 4 Mb with more than 'hotspot_cnt_thresh' number
    # of QTLs)
    
    # Thing to note: Some tests may have chromosomes with multiple hotspots. This results in a problem.
    # One possible solution: only take the bigger of the two for now. 
    # (or largest out of all of them if there are multiple hotspots per chromosomes)
    hotspots = qtl_density_result %>%
      group_by(qtl_chr) %>%
      filter(cnt == max(cnt)) %>%
      filter(cnt > hotspot_cnt_thresh) %>%
      distinct()
    
    # If there are two hotspots taht share a maximum count, get rid of one of them.
    hotspots <- hotspots[!duplicated(hotspots$qtl_chr),]
    hotspots = hotspots %>%
      summarize(center = median(mid)) %>%
      mutate(proximal = center - 2, distal = center + 2)
    
    if (nrow(hotspots) == 0) {next}
    
    # Given the hotspot locations, retain all genes with LOD > threshold
    # and QTL within +/- 4Mb of the mid-point of the hotspot.
    # Returns a list of dataframes, one for each hotspot.
    # Rows of the dataframe contain gene information.
    print("Got hotspot genes...")
    hotspots_genes <- get_hotspots_genes_func(peaks_table, hotspots, qtl_type, lod_threshold, test_type) 
    print("Got hotspot genes...")
    
    # Number of genes in each hotspot.
    hotspots = data.frame(hotspots, count = sapply(hotspots_genes, nrow))
    # kable(hotspots, caption = "Number of genes per hotspot")
    
    ###########################################################################
    ###########################################################################
    # Make plots of correlations between genes in hotspots.
    print("Making hotspot correlation heatmaps...")
    # 1. Plot correlations between their transcipts.
    hotspot_correlations(hotspots_genes, qtl_type = qtl_type, 
                         expr_to_plot = "mrna", test_type = test_type)
    # 2. Plot those gene's corresponding proteins if they have them.
    hotspot_correlations(hotspots_genes, qtl_type = qtl_type, 
                         expr_to_plot = "protein", test_type = test_type)
    
    
    ###########################################################################
    ###########################################################################
    # Make Principal components
    print("Performing principal component analysis...")
    
    hotspots_pcs = as.list(names(hotspots_genes))
    names(hotspots_pcs) = names(hotspots_genes)
    
    pdf(paste0(test_type, "_", qtl_type, "_hotspot_pcs.pdf"), width = 10, height = 8)
    for(i in 1:length(hotspots_genes)) {
      
      # Make title
      chr <- names(hotspots_genes)[i]
      if (grepl("mrna", test_type)) {
        plot_title <- paste0(qtl_type, " eQTLs Chr ", chr, " Gene Hotspot Principal Components (Interactive Covariate = ", covar_name, ")")
      } else {
        plot_title <- paste0(qtl_type, " pQTLs Chr ", chr, " Gene Hotspot Principal Components (Interactive Covariate = ", covar_name, ")")
      }
      
      tmp = hotspots_genes[[i]] %>%
        select(starts_with("DO")) %>%
        as.matrix() %>%
        t() %>%
        pca(method = "svdImpute", nPcs = 3)
      
      hotspots_pcs[[i]] = scores(tmp)
      tmp = gather(data.frame(Mouse.ID = rownames(hotspots_pcs[[i]]), hotspots_pcs[[i]]), key = pc, value = value, -Mouse.ID)
      tmp = left_join(tmp, annot.sample %>% select(Mouse.ID, Age, Sex), by = "Mouse.ID")
      print(tmp %>%
              filter(pc %in% paste0("PC", 1:4)) %>%
              mutate(Age = factor(Age)) %>%
              ggplot(aes(x = Age, value, fill = Sex)) +
              geom_boxplot() +
              facet_grid(pc~.) +
              theme_linedraw() +
              ggtitle(plot_title))
    }
    dev.off()
    
    ###########################################################################
    ###########################################################################
    # Get enriched gene categories for each hotspot. First get the data we are using.
    
    # Pull out data for current expression type.
    # Get type of data. Store column symbols in last row.
    if (grepl("mrna", test_type)) { 
      expr.data <- expr.mrna
      expr.data <- rbind(expr.data, annot.mrna$symbol)
      annot.data <- annot.mrna 
    } else { 
      expr.data <- expr.protein
      expr.data <- rbind(expr.data, annot.protein$symbol)
      annot.data <- annot.protein
    }
    
    # Remove duplicate columns
    expr.data <- expr.data[, !duplicated(colnames(expr.data))]
    
    
    # ###########################################
    # # Allez GSEA
    # ###########################################
    # print("Performing Allez GSEA...")
    # for (i in 1:length(hotspots_genes)) {
    # 
    #   # Subset expr data to only contain expression for genes in our hotspot.
    #   filtered_expr_data <- expr.data[,colnames(expr.data) %in% hotspots_genes[[i]]$lodcolumn]
    #   filtered_expr_data <- as.matrix(filtered_expr_data)
    #   if (ncol(filtered_expr_data) == 1) {next} # Don't bother making a seperate pipleine for hotspots with only one gene. Not worth.
    # 
    #   # Set colnames as gene symbols (gene symbols were stored in the last row) and then remove gene symbols
    #   colnames(filtered_expr_data) <- filtered_expr_data[nrow(filtered_expr_data),]
    #   filtered_expr_data <- filtered_expr_data[-190,]
    # 
    #   # Some gene symbols map to the same ensemble gene id. Thus we might
    #   # have some duplicate gene symbols. Delete duplicates for now
    #   filtered_expr_data <- filtered_expr_data [, !duplicated(colnames(filtered_expr_data ))]
    # 
    # 
    #   ##############
    #   ##############
    #   # Create vector of t statistics for differential expression with respect to either age or sex (depends on test type)
    #   if (grepl("sex", test_type)) {covar = annot.sample$Sex} else {covar = annot.sample$Age}
    # 
    #   score <- vector(length = ncol(filtered_expr_data))
    #   names(score) <- colnames(filtered_expr_data)
    #   for (j in 1:ncol(filtered_expr_data)) {
    #     result <- lm(filtered_expr_data[,j] ~ covar)
    #     t_val <- tidy(result)[2,4] # Pull out t-statistic
    #     score[j] <- t_val
    #   }
    # 
    #   # Perform enrichment analysis. Save results.
    #   result_KEGG <- allez(scores = score, lib = "org.Mm.eg", sets="KEGG", idtype = "SYMBOL")
    #   result_GO <- allez(scores = score, lib = "org.Mm.eg", sets="GO", idtype = "SYMBOL")
    # 
    #   result_KEGG_tbl <- allezTable(result_KEGG)
    #   if (nrow(result_KEGG_tbl) != 0)
    #   {
    #     result_KEGG_tbl <- allez_helper_func(result_KEGG_tbl)
    #     write.csv(result_KEGG_tbl, file = paste0(output_dir, test_type, "_", "chr", names(hotspots_genes)[i],"_", test_type, "_hotspot_KEGG_gsea_allez.csv"))
    #   }
    # 
    #   result_GO_tbl <- allezTable(result_GO)
    #   if (nrow(result_GO_tbl) != 0) {
    #     result_GO_tbl <- allez_helper_func(result_GO_tbl)
    #     write.csv(result_GO_tbl, file = paste0(output_dir, test_type, "_", "chr", names(hotspots_genes)[i],"_", test_type, "_hotspot_GO_gsea_allez.csv"))
    #   }
    # 
    #   # if (nrow(result_GO$setscores) == 0 && nrow(result_KEGG$setscores) == 0) # If both have nothing in them, skip
    #   # {next}
    #   # else if (nrow(result_GO$setscores) == 0) # If only GO has nothing in it, just use KEGG
    #   # {
    #   #   result <- result_KEGG
    #   # }
    #   # else if (nrow(result_KEGG$setscores) == 0) # If only KEGG has nothing in it, just use GO
    #   # {
    #   #   result <- result_GO
    #   # } else { # If both are non-empty, concatenate results.
    #   #   result <- allezC(result_GO, result_KEGG)
    #   # }
    #   # result_tbl <- allezTable(result) # Get significant categories
    #   # if (nrow(result_tbl) == 0) {print("Allez result table has nothing in it!"); next} # If no significant categories, skip.
    # 
    # 
    # } #for (i in 1:length(hotspot_genes)) # Allez

    ###########################################
    # clusterProfiler GSEA
    ###########################################
    print("Performing clusterProfiler GSEA...")
    for (i in 1:length(hotspots_genes)) {
      
      # Subset expr data to only contain expression for genes in our hotspot.
      filtered_expr_data <- expr.data[,colnames(expr.data) %in% hotspots_genes[[i]]$lodcolumn]
      filtered_expr_data <- as.matrix(filtered_expr_data)
      if (ncol(filtered_expr_data) == 1) {next} # Don't bother making a seperate pipleine for hotspots with only one gene. Not worth.
      
      # Set colnames as gene symbols (gene symbols were stored in the last row) and then remove gene symbols
      colnames(filtered_expr_data) <- filtered_expr_data[nrow(filtered_expr_data),]
      filtered_expr_data <- filtered_expr_data[-190,]
      
      # Some gene symbols map to the same ensemble gene id. Thus we might
      # have some duplicate gene symbols. Delete duplicates for now
      filtered_expr_data <- filtered_expr_data [, !duplicated(colnames(filtered_expr_data ))]
      
      
      ##############
      ##############
      # GO over-representation test needs ENTREZID as input
      
      # Create df that maps ENSEMBL, ENTREZ, and symbols to each other.
      # Contains all genes, not just significant genes. (row 190 of expr data contain gene names.)
      universe_df <- bitr(as.character(expr.data[190,]),
                          fromType = "SYMBOL",
                          toType = c("ENTREZID", "ENSEMBL"),
                          OrgDb = org.Mm.eg.db)
      
      
      # Create df that maps ENSEMBL, ENTREZ, and symbols to each other.
      sig_gene_df <- bitr(as.character(colnames(filtered_expr_data)),
                          fromType = "SYMBOL",
                          toType = c("ENTREZID", "ENSEMBL"),
                          OrgDb = org.Mm.eg.db)
      
      
      # ###########################################
      # clusterProfiler Part 1 - GO
      # ###########################################
      # 
      # # GO over-representation test
      go_enrich_output <- enrichGO( gene = sig_gene_df$ENTREZID,
                                    universe = universe_df$ENTREZID,
                                    OrgDb = org.Mm.eg.db,
                                    ont = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05,
                                    readable = TRUE)
      
      go_result <- go_enrich_output@result
      write.csv(go_enrich_output, file = paste0(output_dir, test_type, "_", qtl_type, "_", "chr", names(hotspots_genes)[i], "_hotspot_go_gsea_cluster.csv"), row.names = FALSE)
      
      # ##########################################
      # clusterProfiler Part 2- KEGG
      # ##########################################
      
      # KEGG over-representation test
      kegg_enrich_output <- enrichKEGG(gene = sig_gene_df$ENTREZID,
                                       organism = "mmu",
                                       pvalueCutoff = 0.05)
      if (nrow(kegg_enrich_output) != 0) {
        # Save enrich KEGG, but first convert all ENTREZID into symbol for readability.
        kegg_result <- kegg_enrich_output@result
        kegg_result$geneID_symbol <- NA
        for(k in 1:nrow(kegg_result)) {
          list <- kegg_result$geneID[k]
          list <- strsplit(list, "/")[[1]]
          list <- as.data.frame(list)
          list$ENSEMBL <- NA
          for(n in 1:nrow(list)){
            list$ENSEMBL[n] <- universe_df[universe_df$ENTREZID %in% list[n,"list"],]$SYMBOL[1]
          }
          kegg_result$geneID_symbol[k] <- paste(list$ENSEMBL, collapse = "/")
        }
        write.csv(kegg_result, file = paste0(output_dir, test_type, "_", qtl_type, "_", "chr", names(hotspots_genes)[i], "_hotspot_kegg_gsea_cluster.csv"), row.names = FALSE)
      }
    } #for (i in 1:length(hotspot_genes))
  } #for (qtl_type in c("cis", "trans")) 
} #for (test_type in c("protein_none", "protein_sex", "protein_age","mrna_none", "mrna_sex", "mrna_age")) { 



