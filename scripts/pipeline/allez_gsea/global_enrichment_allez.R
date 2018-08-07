# Date: July 6, 2018
# Author: Daniel Alfonsetti
# Description: Enrichment analysis for DO cross-sectional heart data
# Reference: https://rdrr.io/github/atbroman/allez/

#######################
# Reset environment
rm(list=ls())

# Using local directories
input <- "~/do_heart/data/cleaned_data.RData"
load(input)
output_dir <- "~/do_heart/results/global_enrichment_allez/"
dir.create(output_dir)
# Using HPC directories
# input <- "/home/c-alfond/do_heart/data/cleaned_data.RData"
# load(input)
# output_dir <- "/home/c-alfond/do_heart/scripts/pipeline/output/"

#######################
library(allez)  # Enrichment analysis library
library("org.Mm.eg.db")  # Mouse data base
library(broom)

# Symbol and ID mappings.
ent2ens <- as.list(org.Mm.egENSEMBL)
ent2s <- as.list(org.Mm.egSYMBOL)
#######################
#######################
#######################
helper_func <- function(allez_table) {
  # Description: Takes a table from the output of the allez function and adds columns that contain gene symbols and ensembl ids 
  # corresponding to the set of entrenz gene ids in a given row (category)
  ####
  
  # Debugging
  # allez_table <- sex_result_tbl
  
  # Convert entrenz gene ids into ensemble gene ids and gene symbols
  allez_table[, c("ensembl_ids", "symbols")] <- NA
  names(allez_table)[names(allez_table) == "genes"] <- "entrenz_ids"
  
  # Iterate over each row (category) and update them (add gene ids and gene symbols)
  for (i in 1:nrow(allez_table)) {
    row <- allez_table[i,]
    row$entrenz_ids <- gsub(";", ",", row$entrenz_ids)
    
    entrenz_ids <- strsplit(row$entrenz_ids, ",")[[1]] 
    print(entrenz_ids)
    #result <- select(org.Mm.eg.db, entrenz_ids, c("SYMBOL", "GENENAME"), "ALIAS")
    
    result <- select(org.Mm.eg.db, 
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
#######################
#######################
#######################

# pdf("gsea_plots.pdf", width = 8, height = 6)

# Repeat script for each expression type (mrna and protein)
for (expr.type in c("mrna", "protein")) {
  
  # expr.type <- "protein"
  # Pull out data for current expression type.
  if (expr.type == "mrna") { expr.data <- expr.mrna; colnames(expr.data) <- annot.mrna$symbol; annot.data <- annot.mrna
  } else { expr.data <- expr.protein; colnames(expr.data) <- annot.protein$symbol; annot.data <- annot.protein}
  
  # Remove duplicate columns
  expr.data <- expr.data[, !duplicated(colnames(expr.data))]
  
  # Create vector of t statistics for differential expression with respect to sex
  score_sex <- vector(length = ncol(expr.data))
  names(score_sex) <- colnames(expr.data)
  for (i in 1:ncol(expr.data)) {
    result <- t.test(expr.data[,i] ~ annot.sample$Sex)
    score_sex[i] <- result$statistic # Pull out t-statistic
  }
  
  # Create vector of t statistics for differential expression with respect to age
  score_age <- vector(length = ncol(expr.data))
  names(score_age) <- colnames(expr.data)
  for (i in 1:ncol(expr.data)) {
    result <- lm(expr.data[,i] ~ annot.sample$Age)
    t_val <-tidy(result)[2,4] # Pull out t-statistic
    score_age[i] <- t_val
  }
  
  ###################################################
  ###################################################
  
  # Perform enrichment analysis for sex.
  sex_result_KEGG <- allez(scores = score_sex, lib = "org.Mm.eg", sets="KEGG", idtype = "SYMBOL")
  sex_result_GO <- allez(scores = score_sex, lib = "org.Mm.eg", sets="GO", idtype = "SYMBOL")
  sex_result <- allezC(sex_result_GO, sex_result_KEGG) # Concatenate results
  #allezPlot(sex_result) # Plotting function doesn't work right...
  
  sex_result_tbl <- allezTable(sex_result)
  sex_result_tbl <- helper_func(sex_result_tbl)
  write.csv(sex_result_tbl, file = paste0(output_dir, expr.type, "_sex_gsea_tbl_allez.csv"))    
    
  ###################################################
  
  # Perform enrichment analysis for age. Save results.
  age_result_KEGG <- allez(scores = score_age, lib = "org.Mm.eg", sets="KEGG", idtype = "SYMBOL")
  age_result_GO <- allez(scores = score_age, lib = "org.Mm.eg", sets="GO", idtype = "SYMBOL")
  age_result <- allezC(age_result_GO, age_result_KEGG) # Concatenate results
  #allezPlot(age_result)
  
  age_result_tbl <- allezTable(age_result)
  age_result_tbl <- helper_func(age_result_tbl)
  write.csv(age_result_tbl, file = paste0(output_dir, expr.type, "_age_gsea_tbl_allez.csv"))    
}
#dev.off()
