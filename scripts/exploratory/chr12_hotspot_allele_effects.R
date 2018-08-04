# Author: Daniel Alfonsetti
# Date: July 26, 2018 
# Description: Make QTL plots and look at allele effects stratified by age for 
# interesting genes in our transband.

#######################
# Reset environment
rm(list=ls())
gc()
# Using local directories
# input <- "~/do_heart/data/cleaned_data_with_genoprobs.RData"
# load(input)
# 
# output_dir <- "~/do_heart/results/hotspots/"
# dir.create(output_dir)
# setwd(output_dir)
# num_cores = 2

# Using HPC directories
input <- "/fastscratch/c-alfond/do_heart/data/cleaned_data_with_genoprobs.RData"
load(input)

output_dir <- "/fastscratch/c-alfond/do_heart/results/hotspots/"
dir.create(output_dir)
setwd(output_dir)
num_cores = 1

#######################
# Load in libraries
library(tidyverse)
library(qtl2)
library(doParallel)

registerDoParallel(cores = num_cores)
#######################
#######################

# Setup covariates
rownames(annot.sample) <- unlist(lapply(rownames(annot.sample), gsub, pattern = "-", replacement = "."))

sex <- as.numeric(annot.sample$Sex) # 1 is female, 2 is male
names(sex) <- rownames(annot.sample)

age <- as.numeric(annot.sample$Age)
names(age) <- rownames(annot.sample)

# Old way
# a_covars <- cbind(sex, age, as.numeric(annot.sample$Generation), as.numeric(annot.sample$Batch), as.numeric(as.factor(annot.sample$TMT)))

# New way
df <- annot.sample %>% mutate(TMT <- as.factor(TMT))
a_covars <- model.matrix(~Sex + Age + Generation + Batch + TMT, data = df)[,-1]
rownames(a_covars) <- rownames(annot.sample)
#######################
#######################

# Gene symbol, gene
# Tnt ("ENSMUSP00000011934")
# Myh6 ("ENSMUSP00000080538")
# Obscn ENSMUSP00000049737
# Plec ENSMUSP00000075772
# Myh7b ENSMUSP00000090672
# Dmd ENSMUSP00000109633
# Ccdc22 ENSMUSP00000033483 is a 'control' for diagnostic purposes. It is not in the hotspot. 
foreach(params = list(c("ENSMUSP00000011934", 12), c("ENSMUSP00000080538", 12), c("ENSMUSP00000049737", 12),
                    c("ENSMUSP00000075772", 12), c("ENSMUSP00000090672", 12), c("ENSMUSP00000109633", 12), 
                    c("ENSMUSP00000033483", 12))) %dopar% {
  # for (params in list(c("ENSMUSP00000011934", 12), c("ENSMUSP00000080538", 12))) {
  # Debugging; params <- c("ENSMUSP00000011934", 12) # Tnt
  # params <- c("ENSMUSP00000080538", 12)
  print(params)
  
  id <- params[1]
  qtl_chr <- params[2]
  annot_row <- annot.protein %>% dplyr::filter(protein_id == id)
  protein_id <- id
  symbol <- annot_row$symbol
  expr.data <-  expr.protein[,id]
  
  pdf(paste0("interaction_allele_effects_", toupper(symbol), ".pdf"))
  
  foreach(intcovar_name = c("age", "sex")) %dopar% {
    # Debugging; intcovar_name <- "age"
    intcovar <- eval(parse(text = intcovar_name))
    
    #### Scans
    color <- c("slateblue", "violetred", "green3")
    scan_output_full <- scan1(genoprobs = genoprobs[,qtl_chr], pheno = expr.data, kinship = kin_mat_list[[chr]], addcovar = a_covars, intcovar = intcovar, cores = num_cores)
    scan_output_add <- scan1(genoprobs = genoprobs[,qtl_chr], pheno = expr.data, kinship = kin_mat_list[[chr]], addcovar = a_covars, cores = num_cores)
    scan_output_effect <- scan_output_full - scan_output_add
    
    # Make QTL plot
    plot(x = scan_output_full, map = gmap[qtl_chr], col = color[1])
    plot(x = scan_output_add, map = gmap[qtl_chr], col = color[2], add = TRUE)
    plot(x = scan_output_effect, map = gmap[qtl_chr], col = color[3], add = TRUE)
    
    max_pos <- max(scan_output_full, gmap[qtl_chr])$pos
    if (max_pos > median(gmap[qtl_chr][[1]])) {pos <- "topleft"} else {pos <- "topright"}
    legend(pos, lwd = 2, col = color, c("Full", "Additative", "Interaction Effect"),  
           bg="gray90", lty=c(1,1,1), text.font = 10)
    title(paste0(toupper(symbol), " (", id, ")", " Expression QTL \n with ", intcovar_name, " as an Interactive Covariate"))
    # We want to see place where the effect of the genotype varies with the covariate, not merely
    # where the effect of the QTL was the same for the different values of the covariate.
    # Note: Additative covariates allow to see where 
    
    ##### Run allele effects. Run for all mice and then run stratified by our interactive covariate
    # (either age or sex)
    
    coef <- scan1blup(genoprobs = genoprobs[,qtl_chr],
                      pheno = expr.data,
                      kinship = kin_mat_list[[chr]],
                      addcovar = a_covars)
    plot_coefCC(x = coef, map = gmap[qtl_chr], main = paste0(protein_id))
    
    # Stratify on sex.
    if (intcovar_name == "sex"){

      # Get mouse ids for each sex group.
      annot.sample.M <- annot.sample %>% filter(annot.sample$Sex == "M")
      ids.M <- annot.sample.M$Mouse.ID
      ids.M <- ids.M[ids.M %in% rownames(kin_mat_list[[qtl_chr]])]
      
      
      annot.sample.F <- annot.sample %>% filter(annot.sample$Sex == "F")
      ids.F <- annot.sample.F$Mouse.ID
      ids.F <- ids.F[ids.F %in% rownames(kin_mat_list[[qtl_chr]])]
      
      #######
      
      title <- paste0(toupper(symbol), " (", protein_id, "): Males")
      coef_M <- scan1blup(genoprobs =  genoprobs[,qtl_chr][ids.M],
                          pheno = expr.data[ids.M],
                          kinship = kin_mat_list[[qtl_chr]][ids.M, ids.M],
                          addcovar = a_covars[ids.M,])
      plot_coefCC(x = coef_M, map = gmap[qtl_chr], main = paste0(protein_id, ": Males"))
      
      
      title <- paste0(toupper(symbol), " (", protein_id, "): Females")
      coef_F <- scan1blup(genoprobs =  genoprobs[,qtl_chr][ids.F],
                           pheno = expr.data[ids.F],
                           kinship = kin_mat_list[[qtl_chr]][ids.F, ids.F],
                           addcovar = a_covars[ids.F,])
      plot_coefCC(x = coef_F, map = gmap[qtl_chr], main = paste0(protein_id, ": Females"))
      
    } else { # stratify by age
      
      # Get mouse ids for each age group.
      annot.sample.6 <- annot.sample %>% filter(annot.sample$Age == 6)
      ids.6 <- annot.sample.6$Mouse.ID
      ids.6 <- ids.6[ids.6 %in% rownames(kin_mat_list[[qtl_chr]])]
      
      
      annot.sample.12 <- annot.sample %>% filter(annot.sample$Age == 12)
      ids.12 <- annot.sample.12$Mouse.ID
      ids.12 <- ids.12[ids.12 %in% rownames(kin_mat_list[[qtl_chr]])]
      
      annot.sample.18 <- annot.sample %>% filter(annot.sample$Age == 18)
      ids.18 <- annot.sample.18$Mouse.ID
      ids.18 <- ids.18[ids.18 %in% rownames(kin_mat_list[[qtl_chr]])]
      
      #######

      title <- paste0(toupper(symbol), " (", protein_id, "): 6 months")
      coef_6 <- scan1blup(genoprobs =  genoprobs[,qtl_chr][ids.6],
                          pheno = expr.data[ids.6],
                          kinship = kin_mat_list[[qtl_chr]][ids.6, ids.6],
                          addcovar = a_covars[ids.6,])
      plot_coefCC(x = coef_6, map = gmap[qtl_chr], main = paste0(protein_id, ": 6 months"))
      
      
      title <- paste0(toupper(symbol), " (", protein_id, "): 12 months")
      coef_12 <- scan1blup(genoprobs =  genoprobs[,qtl_chr][ids.12],
                           pheno = expr.data[ids.12],
                           kinship = kin_mat_list[[qtl_chr]][ids.12, ids.12],
                           addcovar = a_covars[ids.12,])
      plot_coefCC(x = coef_12, map = gmap[qtl_chr], main = paste0(protein_id, ": 12 months"))
      
      
      title <- paste0(toupper(symbol), " (", protein_id, "): 18 months")
      coef_18 <- scan1blup(genoprobs =  genoprobs[,qtl_chr][ids.18],
                           pheno = expr.data[ids.18],
                           kinship = kin_mat_list[[qtl_chr]][ids.18, ids.18],
                           addcovar = a_covars[ids.18,]) 
      plot_coefCC(x = coef_18, map = gmap[qtl_chr], main = title)
      
    } # if (intcovar_name == "age"){
  } # (intcovar_name in c("age", "sex")) {
  dev.off()
} # (params in c( ...