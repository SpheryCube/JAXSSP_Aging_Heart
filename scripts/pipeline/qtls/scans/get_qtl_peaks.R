################################################################################
# Generate eQTL and pQTL peak tables with allele effects for each test type.
# R/3.4.4
# Daniel Alfonsetti
# daniel.alfonsetti@gmail.com
# Jul. 22, 2018
################################################################################
# Reset environment
rm(list = ls())

library(qtl2)
library(tidyverse)
library(doParallel)
library(foreach)

######## 
# Using Local directories
######## 
# num_cores <- 1
# lod_thresh <- 6
# 
# # Input paths
# genoprobs_file <- "~/do_heart/data/cleaned_data_with_genoprobs.RData"
# load(genoprobs_file)
# 
# perms_input_dir <- "~/do_heart/results/"
# qtl_input_dir <- "~/do_heart/results/"
# 
# # # Output paths
# peaks_output_dir <- "~/do_heart/results/"

######## 
# Using HPC directories
######## 
num_cores <- 12
lod_thresh <- 6    # Keep this at 6, but set threshold for QTL hotspots at 7.2

# Input paths
genoprobs_file <- "/home/c-alfond/do_heart/data/cleaned_data_with_genoprobs.RData"
load(genoprobs_file)

# perms_input_dir <- "/fastscratch/c-alfond/do_heart/scripts/pipeline/output/" # In case we are using perms threshold.
qtl_input_dir <- "/fastscratch/c-alfond/do_heart/results/"

# Output paths
peaks_output_dir <- "/fastscratch/c-alfond/do_heart/results/"
##########################################################
##########################################################

# Setup some varables needed for the coef_fxn
sex <- as.numeric(annot.sample$Sex) # 1 is female, 2 is male
names(sex) <- rownames(annot.sample)

age <- as.numeric(annot.sample$Age)
names(age) <- rownames(annot.sample)

# Register parallele back end.
registerDoParallel(cores = num_cores)
#######################
#######################

# Helper function
# For each peak, get the BLUP coefficients and add them to the output
coef_fxn = function(peaks, expr_type) {
  
  # Setup variables according to parameters
  # Include sex, age, generation, batch, and tag for protein expression data. 
  # Include only sex, age, and generation for mrna.
  if (expr_type == 'protein')
  {
    pheno_data <- expr.protein
    
    # Old
    # a_covars <- cbind(sex, age, as.numeric(annot.sample$Generation), as.numeric(annot.sample$Batch), as.numeric(as.factor(annot.sample$TMT)))
    # New
    # df <- annot.sample %>% mutate(TMT <- as.factor(TMT))
    # a_covars <- model.matrix(~Sex + Age + Generation + Batch + TMT, data = df)[,-1]
    # rownames(a_covars) <- rownames(annot.sample)
    # New 2. Don't want to overfit.
    a_covars <- cbind(sex, age)

   } else {
    pheno_data <- expr.mrna
    
    # Old
    # a_covars <- cbind(sex, age, as.numeric(annot.sample$Generation))
    # New
    # a_covars <- model.matrix(~Sex + Age + Generation, data = annot.sample)[,-1]
    # rownames(a_covars) <- rownames(annot.sample)
    # New 2. Don't want to overfit.
    a_covars <- cbind(sex, age)
  }
  
  # Add columns for allele effect coefficients
  peaks = cbind(peaks, matrix(0, nrow = nrow(peaks), ncol = 8, dimnames =
                                list(NULL, LETTERS[1:8])))
  cols2fill = (ncol(peaks) - 7):ncol(peaks)
  
  # Iterate over rows in parallel
  peaks_w_effects <- foreach(i = 1:nrow(peaks), .combine = "rbind") %dopar% {
    qtl_chr  = peaks$qtl_chr[i]
    qtl_pos <- peaks$qtl_pos[i]
    
    # Get marker id
    chr_gmap <- gmap[[qtl_chr]]
    tmp <- which.min(abs(chr_gmap-qtl_pos))
    qtl_mkr <- names(tmp)
    rm(tmp)

    # ensembl id for expression data being used. Could be protein or gene id.
    pheno_id = peaks$lodcolumn[i] 
    
    gp = genoprobs[,qtl_chr] # Get genoprobs for this chromosome at this marker.
    gp[[1]] = gp[[1]][,,qtl_mkr, drop = FALSE] 
    
    # Additive scans
    blup = scan1blup(genoprobs = gp, pheno = pheno_data[,pheno_id, drop = FALSE],
                     kinship = kin_mat_list[[qtl_chr]], addcovar = a_covars, 
                     cores = num_cores)

    peaks[i,cols2fill] = blup[1,1:8]
    
    peaks[i,] # return row
  } # foreach()
  return(peaks_w_effects)
} # coef_fxn()

#######################
#######################

# Extract peak tables. Do in parallel
foreach(scan_type = c("protein_age", "protein_sex", "protein_none", "protein_age_full", "protein_sex_full")) %dopar% {
# foreach (scan_type = c("mrna_age_full", "mrna_sex_full")) %dopar% {
  print("=====================================")
  print(scan_type)

  # Read in scan
  print("Loading in scan...")
  qtl_file_name <- paste0("qtlscan_", scan_type, ".rds")
  qtl_scan <- readRDS(paste0(qtl_input_dir, qtl_file_name))

  print("Getting peaks!")

  # If we are using perms....
    # Get current scans' corresponding set of permutations.
    # perm_file_name <- paste0("perms_", scan_type, ".rds")
    # perms <- readRDS(paste0(perms_input_dir, perm_file_name))
    #
    # # # Get 95 percentile of the permutation LOD distribution for each phenotype(expression data for an mRNA transcript or protein),
    # # using the permutation object.
    # thresholds95 <- summary(perms)
    # peaks_table <- find_peaks(scan1_output = qtl_scan, map = gmap, threshold = thresholds95, prob = 0.95, cores = num_cores)

  # If we are using a flat cutoff...
  peaks_table <- find_peaks(scan1_output = qtl_scan, map = gmap, threshold = lod_thresh, prob = 0.95, cores = num_cores)
  
  if (nrow(peaks_table) == 0) {print(paste0("0 peaks for: ", scan_type)); next}
  
  print("Got peaks!")
  # Rename columns
  peaks_table <- dplyr::rename(peaks_table, qtl_chr = chr, qtl_pos = pos, qtl_lod = lod)

  # Add columns
  namevector <- c("gene_symbol","gene_chr","gene_start", "gene_middle", "gene_end", "cis", "qtl_mkr")
  peaks_table[, namevector] <- NA

  levels <- c(1:19, "X", "Y", "MT")
  peaks_table$gene_chr <- factor(peaks_table$gene_chr, levels)

  # Add gene symbol and gene chromosome for each QTL (and gene id for pQTLs specifically)
  print("Adding gene symbol and gene chromosome for each QTL...")
  if (grepl("mrna", scan_type)) { # we are dealing with an mrna scan

    # Reoder columns
    peaks_table <- peaks_table[, c("lodcolumn", "gene_symbol","gene_chr", "gene_start", "gene_middle",
                                   "gene_end", "qtl_chr", "qtl_pos", "qtl_mkr", "qtl_lod", "cis")]

    # Update each row
    peaks_table <- foreach(i = 1:nrow(peaks_table), .combine = "rbind") %dopar% {
      row <- peaks_table[i,]

      annot.row <- annot.mrna %>% filter(gene_id == row$lodcolumn)
      annot.row <- annot.row[1,]  # Make sure we are only pulling out one row. (some transcripts map to multiple proteins)

      row$gene_symbol <- annot.row$symbol
      row$gene_chr <- annot.row$chr
      row$gene_start <- annot.row$start / 10^6 # Convert bases to mega bases
      row$gene_middle <- annot.row$middle / 10^6
      row$gene_end <- annot.row$end / 10^6

      peaks_table[i,] <- row
      peaks_table[i,]
    }

  } else {  # Basically the same logic, except we are dealing with a protein scan.

    # Add an extra column for gene_id (since protein ids are already in the "lodcolumn" column.)
    peaks_table[,"gene_id"] <- NA
    # Reorder columns
    peaks_table <- peaks_table[, c("lodcolumn", "gene_id", "gene_symbol","gene_chr", "gene_start", "gene_middle",
                                   "gene_end", "qtl_chr", "qtl_pos", "qtl_mkr", "qtl_lod", "cis")]

    # Update each row
    peaks_table <- foreach(i = 1:nrow(peaks_table), .combine = "rbind") %dopar% {
      row <- peaks_table[i,]

      annot.row <- annot.protein %>% filter(protein_id == row$lodcolumn)
      annot.row <- annot.row[1,]

      row$gene_symbol <- annot.row$symbol
      row$gene_chr <- annot.row$chr
      row$gene_start <- annot.row$start / 10^6
      row$gene_middle <- annot.row$middle / 10^6
      row$gene_end <- annot.row$end / 10^6

      row$gene_id <- annot.row$gene_id

      peaks_table[i,] <- row
      peaks_table[i,]
    }
  }

  # Sort by LOD score.
  print("Soring by LOD score...")
  peaks_table <- peaks_table[order(peaks_table$qtl_lod, decreasing =TRUE),]

  # Add allele effects
  # Can only get allele effects for models that only use additive covariates.
  if (scan_type == "protein_none" | scan_type == "mrna_none") {
    print("Getting allele effects...")
    expr_type <- ifelse(grepl("mrna", scan_type), "mrna", "protein")
    peaks_table <- coef_fxn(peaks_table, expr_type)
    print("Got allele effects!")
  }

  # Add QTL marker ids
  print("Adding QTL marker ids")
  peaks_table <- foreach(i = 1:nrow(peaks_table), .combine = "rbind") %dopar% {
    qtl_chr  = peaks_table$qtl_chr[i]
    qtl_pos <- peaks_table$qtl_pos[i]

    chr_gmap <- gmap[[qtl_chr]]
    tmp <- which.min(abs(chr_gmap-qtl_pos))
    qtl_mkr <- names(tmp)
    peaks_table[i, "qtl_mkr"] <- qtl_mkr
    rm(tmp)
    peaks_table[i,] # return row
  }

  # Set levels for qtl_chr and gene_chr. Force into factors.
  print("Setting factor levels for qtl_chr and gene_chr...")
  levels <- c(1:19, "X", "Y", "MT")
  peaks_table$qtl_chr <- factor(peaks_table$qtl_chr, levels = levels)
  peaks_table$gene_chr <- factor(peaks_table$gene_chr, levels = levels)

  # Display whether the qtl is local or distal.
  # QTLs are considered local ('cis') if they are within 2 megabases of each other on the same chromosome.
  peaks_table[, 'cis'] <- (peaks_table$gene_chr == peaks_table$qtl_chr) & (abs(peaks_table$gene_middle - peaks_table$qtl_pos) < 4)

  # Save peaks table
  print("Saving...")
  peaks_file_name <- paste0("peaks_", scan_type)
  write.csv(peaks_table, file = paste0(peaks_output_dir, peaks_file_name, "_thresh_", as.character(lod_thresh), ".csv"), row.names = FALSE)
  print("Done making allele tables")
}

####################################################################################
####################################################################################
####################################################################################
# Diagnostics/sanity checks for the mrna_none and protein_none tests
# Test a few sets of coefficients from our outputs to make sure
# they are congruent with effect plots.

# setwd(peaks_output_dir)
# print("Starting diagnostics...")
# print("Diagnostics for protein_none...")
# pdf("peaks_diagnostics.pdf")
# 
# # Test first for protein expression level allele effects.
# scan_type <- "protein_none"
# 
# qtl_file_name <- paste0("qtlscan_", scan_type, ".rds")
# qtl_scan <- readRDS(paste0(qtl_input_dir, qtl_file_name))
# a_covars <- cbind(sex, age, as.numeric(annot.sample$Generation), as.numeric(annot.sample$Batch), as.numeric(as.factor(annot.sample$TMT)))
# rownames(a_covars) <- unlist(lapply(rownames(a_covars), gsub, pattern = "-", replacement = "."))
# 
# # ENSMUSP00000139261 (KCTD12): QTL on chr 4 at 105 Mb. 
# protein_id <- "ENSMUSP00000139261"
# chr <- 4
# 
# coef <- scan1blup(genoprobs = genoprobs[,chr], 
#                   pheno = expr.protein[,protein_id], 
#                   kinship = kin_mat_list[[chr]],
#                   addcovar = a_covars)
# scan_output <- scan1(genoprobs = genoprobs, pheno = expr.protein[,protein_id], kinship = kin_mat_list, addcovar = a_covars, cores = num_cores)
# plot_coefCC(x = coef, map = gmap[chr], scan1_output = scan_output, main = protein_id)
# 
# # ENSMUSP00000138597: QTl on chr 7 at 134 Mb. Gga3
# protein_id <- "ENSMUSP00000138597"
# chr <- 7
# 
# coef <- scan1blup(genoprobs = genoprobs[,chr], 
#                   pheno = expr.protein[,protein_id], 
#                   kinship = kin_mat_list[[chr]],
#                   addcovar = a_covars)
# scan_output <- scan1(genoprobs = genoprobs, pheno = expr.protein[,protein_id], kinship = kin_mat_list, addcovar = a_covars, cores = num_cores)
# plot_coefCC(x = coef, map = gmap[chr], scan1_output = scan_output, main = protein_id)
# 
# # ENSMUSP00000096753: Nnt qtl_chr 13 at qtl_pos 118
# protein_id <- "ENSMUSP00000096753"
# chr <- 13
# 
# coef <- scan1blup(genoprobs = genoprobs[,chr], 
#                   pheno = expr.protein[,protein_id], 
#                   kinship = kin_mat_list[[chr]],
#                   addcovar = a_covars)
# scan_output <- scan1(genoprobs = genoprobs, pheno = expr.protein[,protein_id], kinship = kin_mat_list, addcovar = a_covars, cores = num_cores)
# plot_coefCC(x = coef, map = gmap[chr], scan1_output = scan_output, main = protein_id)
# 
# #############
# # Now do diagnostics for mrna
# ############
# print("Diagnostics for mrna_none...")
# 
# scan_type <- "mrna_none"
# 
# qtl_file_name <- paste0("qtlscan_", scan_type, ".rds")
# qtl_scan <- readRDS(paste0(qtl_input_dir, qtl_file_name))
# a_covars <- cbind(sex, age, as.numeric(annot.sample$Generation))
# rownames(a_covars) <- unlist(lapply(rownames(a_covars), gsub, pattern = "-", replacement = "."))
#  
# 
# #ENSMUSG00000099162 (RP23-279M5.2). QTL on chr X at 131 Mb
# gene_id <- "ENSMUSG00000099162"
# chr <- "X"
# 
# coef <- scan1blup(genoprobs = genoprobs[,chr], 
#                   pheno = expr.mrna[,gene_id], 
#                   kinship = kin_mat_list[[chr]],
#                   addcovar = a_covars)
# scan_output <- scan1(genoprobs = genoprobs, pheno = expr.mrna[,gene_id], kinship = kin_mat_list, addcovar = a_covars, cores = num_cores)
# plot_coefCC(x = coef, map = gmap[chr], scan1_output = scan_output, main = gene_id)
# 
# # 
# # ENSMUSG00000041453	Rpl21, QTL on chr 5 at 146. QTL lod of 83
# gene_id <- "ENSMUSG00000041453"
# chr <- 5
# 
# coef <- scan1blup(genoprobs = genoprobs[,chr], 
#                   pheno = expr.mrna[,gene_id], 
#                   kinship = kin_mat_list[[chr]],
#                   addcovar = a_covars)
# scan_output <- scan1(genoprobs = genoprobs, pheno = expr.mrna[,gene_id], kinship = kin_mat_list, addcovar = a_covars, cores = num_cores)
# plot_coefCC(x = coef, map = gmap[chr], scan1_output = scan_output, main = gene_id)
# dev.off()
# print("Done.")
# 
# # Diagnostics look good.