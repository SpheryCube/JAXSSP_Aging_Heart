# Author: Daniel Alfonsetti (daniel.alfonsetti@gmail.com)
# Date: July 31, 2018

# For each hotspot, look for mediators underneath.
# Since we are doing interactive QTLs, we look to see if the mediator gets rid of the interaction effect.

# [y ~ intcovar + Q + Q:intcovar] -  [y ~ intcovar + Q] = interaction effect
# [y ~ M + intcovar + Q + Q:intcovar] - [y ~ M + intcovar + Q] = interaction effect with candidate mediator as covariate
# If the interaction effect with the candidate mediator is much smaller than the interaction effect w/o the candidate mediator,
# then we know the candidate mediator is in fact a mediator and is mediating the interactive effect.

#######################
rm(list = ls())
gc()
library(qtl2)
library(tidyverse)
library(doParallel)

# Local directories
# load("~/do_heart/data/cleaned_data_with_genoprobs.RData")
# query_genes = create_gene_query_func(dbfile = "~/do_heart/data/mouse_genes.sqlite", filter = "source='MGI'")
# output_dir <- "~/do_heart/results/hotspots"
# setwd(output_dir)
# base_hotspots_dir <- "/Users/c-alfond/do_heart/results/hotspots/"
# 
# 
# qtl_dir <- "~/do_heart/results/"
# num_cores <- 1

#######################
# HPC Directories
load("/fastscratch/c-alfond/do_heart/data/cleaned_data_with_genoprobs.RData")
query_genes = create_gene_query_func(dbfile = "/fastscratch/c-alfond/do_heart/data/mouse_genes.sqlite", filter = "source='MGI'")
output_dir <- "/fastscratch/c-alfond/do_heart/results/"
setwd(output_dir)
base_hotspots_dir <- "/fastscratch/c-alfond/do_heart/results/hotspots/"

num_cores <- 35
radius <- 10

# #################### 
# Register parallele back end.
registerDoParallel(cores = num_cores)
# #################### 
# ####################
sex <- as.numeric(annot.sample$Sex) # 1 is female, 2 is male
names(sex) <- rownames(annot.sample)

age <- as.numeric(annot.sample$Age)
names(age) <- rownames(annot.sample)

a_covars <- cbind(sex, age)

#######################


for (test_type in c("mrna_sex")) {
# foreach (test_type = c("protein_age", "protein_sex", "mrna_age", "mrna_sex")) %dopar% {
  # foreach (test_type = c("mrna_sex")) %dopar% {
  # Debugging
    # test_type <- "protein_age"
  #
  print(test_type)
    
  if (grepl("protein", test_type)){
    intcovar <- age
  } else {
    intcovar <- sex
  }

  
  foreach (mediation_type = c("mrna", "protein")) %dopar% {
  # foreach (mediation_type = c("protein")) %dopar% {
  # for (mediation_type in c("mrna", "protein")) {
    # Debugging
      # mediation_type <- "protein"
    #
    print(paste0("Mediation_type: ", mediation_type))
    
    # Intialize a dataframe to store result for genes that are mediators
    
    mediators_df <- data.frame(matrix(ncol=15, nrow = 0))
    names(mediators_df) <- c("target_id", "target_symbol", "gene_chr", "gene_pos", "qtl_chr", "qtl_pos", "mediator_id", "mediator_symbol",
                                         "mediator_chr", "mediator_pos", "intscan_w_candidate", "addscan_w_candidate", "int_effect_w_candidate", "int_effect_lod", "diff")
    
    
    setwd(paste0(base_hotspots_dir, test_type, "_hotspots_output/"))
    hotspot_gene_files <- dir(pattern ="*hotspot_genes.csv")
    
    for (file in hotspot_gene_files) {
      # Debugging
        # file <- "protein_age_trans_chr3_hotspot_genes.csv"
      #
      hotspot_genes_df <- read.csv(file)
      
      # Debugging. Saves time
        # hotspot_genes_df <- hotspot_genes_df[1:2, ]
      for (i in 1:nrow(hotspot_genes_df)) {

        target_gene_row <- hotspot_genes_df[i,]
        target_symbol <- as.character(target_gene_row$gene_symbol)
        target_id <- as.character(target_gene_row$lodcolumn)
        target_gene_chr <- target_gene_row$gene_chr
        target_gene_pos <- target_gene_row$gene_middle
        
        qtl_chr <- target_gene_row$qtl_chr
        qtl_pos <- target_gene_row$qtl_pos
        int_effect_lod <- target_gene_row$qtl_lod 
        max_qtl_mkr <- as.character(target_gene_row$qtl_mkr)
        
        # Interval to search for mediators in.
        start <- qtl_pos - radius
        end <- qtl_pos + radius
        
        # Get the expression type of the target gene.
        if (grepl("protein", test_type)) {
          data.target <- expr.protein[,target_id]
        } else {
          data.target <- expr.mrna[,target_id]
        }
        
        # Get candidate mediators for target gene These will be genes under the hotspot.
        genes = query_genes(qtl_chr, start, end)
        candidate_gene_names <- genes$Name
        str(genes)
        if (nrow(genes) == 0){next}
        
        if (grepl("protein", mediation_type)){
          candidate_genes_annots <- annot.protein[annot.protein$symbol %in% candidate_gene_names,, drop = FALSE]
        } else {
          candidate_genes_annots <- annot.mrna[annot.mrna$symbol %in% candidate_gene_names,, drop = FALSE]
        }
        if (nrow(candidate_genes_annots) == 0){next}
        

        for (k in 1:nrow(candidate_genes_annots)) {
          candidate_gene_row <- candidate_genes_annots[k,]
          
          # Get expression date depending on what type (mRNA or protein) we are mediating on.
          if ("protein" == mediation_type) {
            mediator_id <- candidate_gene_row$protein_id
            candidate_gene_expr <- expr.protein[, mediator_id, drop = FALSE]
          } else {
            mediator_id <- candidate_gene_row$gene_id
            candidate_gene_expr <- expr.mrna[, mediator_id, drop = FALSE]
          }
          mediator_symbol <- candidate_gene_row$symbol
          mediator_chr <- candidate_gene_row$chr
          mediator_pos <- candidate_gene_row$middle / 10^6
          
          a_covars_filtered <- a_covars[rownames(a_covars) %in% rownames(candidate_gene_expr),]
          a_covars_filtered_w_candidate <- merge(a_covars_filtered, candidate_gene_expr, by = "row.names")
          rownames(a_covars_filtered_w_candidate) <- a_covars_filtered_w_candidate[,1]
          a_covars_filtered_w_candidate <- a_covars_filtered_w_candidate[,-1]
          a_covars_filtered_w_candidate <- a_covars_filtered_w_candidate[complete.cases(a_covars_filtered_w_candidate),]
          
          
          ##########################
          ##########################
          ##########################
          # intscan_wo_candidate <- scan1(genoprobs  = genoprobs[,qtl_chr], pheno = data.target, kinship = kin_mat_list[[qtl_chr]],
          #                               intcovar = intcovar, addcovar = a_covars)
          # intscan_wo_candidate[max_qtl_mkr,1]
          # 
          # addscan_wo_candidate <-  scan1(genoprobs  = genoprobs[,qtl_chr], pheno = data.target, kinship = kin_mat_list[[qtl_chr]],
          #                               addcovar = a_covars)
          # addscan_wo_candidate[max_qtl_mkr,1]
          
          ##########################
          ##########################
          intscan_w_candidate <- scan1(genoprobs  = genoprobs[,qtl_chr], pheno = data.target, kinship = kin_mat_list[[qtl_chr]], 
                                       intcovar = intcovar, addcovar = a_covars_filtered_w_candidate, num_cores = 15)
          intscan_w_candidate[max_qtl_mkr,1]
          

          addscan_w_candidate <-  scan1(genoprobs  = genoprobs[,qtl_chr], pheno = data.target, kinship = kin_mat_list[[qtl_chr]], 
                                        addcovar = a_covars_filtered_w_candidate, num_cores = 15)
          addscan_w_candidate[max_qtl_mkr,1]
          
          
          # Get effect of intersection with models that have candidate mediator
          int_effect_w_candidate <- intscan_w_candidate[max_qtl_mkr,1] -  addscan_w_candidate[max_qtl_mkr,1]
          
          # Compare effects. (Essentially a difference of differences)
          diff <- int_effect_lod - int_effect_w_candidate
          
          # Round  
          int_effect_w_candidate <- round(int_effect_w_candidate, 4)
          int_effect_lod <- round(int_effect_lod, 4)
          diff <- round(diff, 4)
          
          # Append mediation results to output df.
          new_row <- c(target_id, target_symbol, target_gene_chr, target_gene_pos, qtl_chr, qtl_pos, mediator_id, 
                       mediator_symbol, mediator_chr, mediator_pos, intscan_w_candidate[max_qtl_mkr,1], addscan_w_candidate[max_qtl_mkr,1], 
                       int_effect_w_candidate, int_effect_lod, diff)
          mediators_df[nrow(mediators_df) + 1,] <- new_row

        } #  for (k in 1:nrow(candidate_genes_annots)) {
      } # for (i in 1:nrow(hotspot_genes_df)) {
    }# for (file in hotspot_gene_files) {
    
    mediators_df <- mediators_df[order(-as.numeric(mediators_df$diff)),]
    write.csv(mediators_df, file = paste0(test_type, "_hotspot_", mediation_type, "_mediators_r", radius,".csv"), row.names = FALSE)
  }   # for (mediation_type in c("mrna", "protein"))
} # for (test_type in c("protein_age", "protein_sex", "mrna_age", "mrna_sex"))

    


