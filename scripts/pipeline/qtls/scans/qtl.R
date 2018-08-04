# Author: Daniel Alfonsetti (daniel.alfonsetti@gmail.com)
# Date: July 20th, 2018
# Description: Run pQTLs and eQTLs in parallel on all protein and mRNA expressions and then again considering age and sex as interactive covariates.

# References: 
# http://kbroman.org/qtl2/assets/vignettes/user_guide.html#qtl_analysis_in_diversity_outbred_mice
# https://smcclatchy.github.io/mapping/13-qtl-in-do/

#######################
# Reset environment
rm(list = ls())

library(qtl2)
library(doParallel)
library(foreach)
library(dplyr)



# Using local directories
# load("~/do_heart/data/cleaned_data_with_genoprobs.RData")
# output_dir <- "~/do_heart/results/"
# 
# num_cores <- 3
# num_perms <- 1


# Using HPC directories
load("/home/c-alfond/do_heart/data/cleaned_data_with_genoprobs.RData")
output_dir <- "/fastscratch/c-alfond/do_heart/results/"
num_cores <- 10
num_perms <- 10000

# #################### 
# Register parallele back end.
registerDoParallel(cores = num_cores)

#########################################################
# Setup covariate objects
sex <- as.numeric(annot.sample$Sex) # 1 is female, 2 is male
names(sex) <- rownames(annot.sample)

age <- as.numeric(annot.sample$Age)
names(age) <- rownames(annot.sample)

none <- NA

####################

# Run all possible combinations of paramters in parallel.
foreach(p_expr_type = c("protein", "mrna")) %dopar% {
  foreach(p_covar = c("age", "sex", "none")) %dopar% {
      print(paste0("Expression type is ", p_expr_type, " and the interactive covariate is ", p_covar))

      # Extract objects according to the parameters for this job.
      expr_data <- eval(parse(text = paste0('expr.', p_expr_type))) # Either mrna or protein
      covar <- eval(parse(text = p_covar)) # Either age, sex, or none. (Determines interactive covariate type)
      
      # Debugging
        # expr_data <- expr_data[,1:2]
      print(expr_data)
      #####################

      # Create additive covariates to reduce noise.
        # This was commented out in fear that we were overfitting the data before.
        # The data is already normalized so we shouldn't ahve to add batch, tag, and generation covars.
          # Include sex, age, generation, batch, and tag for protein. Include only sex, age, and generation for mrna.
          # if (p_expr_type == 'protein')
          # {
          #   # a_covars <- cbind(sex, age, as.numeric(annot.sample$Generation), as.numeric(annot.sample$Batch), as.numeric(as.factor(annot.sample$TMT)))
          #   
          #   df <- annot.sample %>% mutate(TMT <- as.factor(TMT))
          #   a_covars <- model.matrix(~Sex + Age + Generation + Batch + TMT, data = df)[,-1]
          #   rownames(a_covars) <- rownames(annot.sample)
          #   
          # } else {
          #   # a_covars <- cbind(sex, age, as.numeric(annot.sample$Generation))
          # 
          #   a_covars <- model.matrix(~Sex + Age + Generation, data = annot.sample)[,-1]
          #   rownames(a_covars) <- rownames(annot.sample)
          # }
      a_covars <- cbind(sex, age)

      # Run scan and corresponding permutation.
      # Because the data is normalized, we only have to do a permutation testing on one gene per scan type and then apply the threshold to the rest of the genes.
      if (p_covar == 'none') {
        print(paste0("Scanning with ", p_expr_type, " ", p_covar))
        scan_output <- scan1(genoprobs = genoprobs, pheno = expr_data, kinship = kin_mat_list, addcovar = a_covars, cores = num_cores)
        # print("Permuting")
        # perms_output <- scan1perm(genoprobs = genoprobs, pheno = expr_data[,1], kinship = kin_mat_list, addcovar = a_covars, cores = num_cores, n_perm = num_perms)
      }

      if (p_covar == 'sex' | p_covar == 'age')
      {
        print(paste0("Scanning with ", p_expr_type, " ", p_covar))
        scan_output <- scan1(genoprobs = genoprobs, pheno = expr_data, kinship = kin_mat_list, addcovar = a_covars, intcovar = covar, cores = num_cores)
        # print("Permuting")
        # perms_output <- scan1perm(genoprobs = genoprobs, pheno = expr_data[,1], kinship = kin_mat_list, addcovar = a_covars, intcovar = covar, cores = 5, n_perm = num_perms)
      }

      # Save results of scan and permutation
      print("Saving full scan...")
      if (p_covar == "none") {
        scan_file_name <- paste0("qtlscan_", p_expr_type, "_", p_covar, ".rds")
      } else {
        scan_file_name <- paste0("qtlscan_", p_expr_type, "_", p_covar, "_full.rds")
      }

      saveRDS(scan_output, file = paste0(output_dir, scan_file_name))
      # print("Saving permutations...")
      # perm_file_name <- paste0("perms_", expr_type, "_", covar, "_full.rds")
      # saveRDS(perms_output, file = paste0(output_dir, perm_file_name))
      
  }
}


# Get interaction effects using qtl scans we got from above.
print("Getting interaction effects...")
setwd(output_dir)
files <- dir(pattern = "qtlscan")

files1 <- files[grep("full", files)]
files2 <- files[grep("none", files)]
files <- c(files1, files2)

foreach(p_expr_type = c("mrna", "protein")) %dopar% {
  # debugging p_expr_type <- "protein"
  expr_files <- files[grepl(p_expr_type, files)]
  print(expr_files)
  expr_none_file <- expr_files[grepl("none", expr_files)]
  print(expr_none_file)
  expr_none <- readRDS(expr_none_file)

  foreach(p_covar = c("age", "sex")) %dopar% {
    # Debugging, p_covar <- "age"
    print(paste0("Getting interaction effects for ", p_expr_type, " ", p_covar))
    expr_p_covar_file <- expr_files[grepl(p_covar, expr_files)]
    expr_p_covar <- readRDS(expr_p_covar_file)
    
    interaction_result <- expr_p_covar - expr_none
    print("Saving interaction scan")
    saveRDS(interaction_result, file = paste0("qtlscan_", p_expr_type, "_", p_covar, "_test.rds"))
  }
}

