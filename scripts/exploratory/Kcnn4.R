# Author: Daniel Alfonsetti
# Date: July 30, 2018 
# We want to see place where the effect of the genotype varies with the covariate

#######################
# Load in libraries
library(tidyverse)
library(qtl2)
library(doParallel)
##############################################
##############################################
# # Using local directories
# rm(list=ls())
# gc()
# input <- "~/do_heart/data/cleaned_data_with_genoprobs.RData"
# load(input)
# 
# output_dir <- "~/do_heart/results/hotspots/"
# hotspots_dir <- "~/do_heart/results/hotspots/"
# mouse_genes.sqlite <- "~/do_heart/data/mouse_genes.sqlite"
# query_func = create_variant_query_func("~/do_heart/data/cc_variants.sqlite")
# 
# num_cores = 2
###################
# Using HPC directories
rm(list=ls())
gc()
input <- "/fastscratch/c-alfond/do_heart/data/cleaned_data_with_genoprobs.RData"
load(input)
#
output_dir <- "/fastscratch/c-alfond/do_heart/results/hotspots/"
hotspots_dir <- "/fastscratch/c-alfond/do_heart/results/hotspots/"
mouse_genes.sqlite <- "/fastscratch/c-alfond/do_heart/data/mouse_genes.sqlite"
query_func = create_variant_query_func("/fastscratch/c-alfond/do_heart/data/cc_variants.sqlite")

num_cores = 38
##############################################
##############################################
registerDoParallel(cores = num_cores)

# Setup covariates
rownames(annot.sample) <- unlist(lapply(rownames(annot.sample), gsub, pattern = "-", replacement = "."))

sex <- as.numeric(annot.sample$Sex) # 1 is female, 2 is male
names(sex) <- rownames(annot.sample)

age <- as.numeric(annot.sample$Age)
names(age) <- rownames(annot.sample)

color <- c("slateblue", "violetred", "green3") # for QTLs

##############################################
##############################################
##############################################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
##############################################
##############################################
expr_v_age_by_genotype_plots <- function(expr.data, genotypes_snp, snpinfo, annot.sample, symbol, id) {
  
  expr.data <- as.data.frame(expr.data)
  
  # Setup covars
  sexF <- annot.sample$Sex # factor
  names(sexF) <- rownames(annot.sample)
  
  age <- annot.sample$Age # numeric
  names(age) <- rownames(annot.sample)
  
  # Append sex, age, and expression data to genotypes.
  data.w.max.snp.genos <- merge(expr.data,  genotypes_snp, by = "row.names")
  rownames(data.w.max.snp.genos) <- data.w.max.snp.genos[,"Row.names"] # Reset rownames.
  data.w.max.snp.genos <- data.w.max.snp.genos[,-1] # Delete dedicated rownames column.
  
  data.w.max.snp.genos <- merge(data.w.max.snp.genos,  sexF, by = "row.names")
  rownames(data.w.max.snp.genos) <- data.w.max.snp.genos[,"Row.names"]
  data.w.max.snp.genos <- data.w.max.snp.genos[,-1]
  
  data.w.max.snp.genos <- merge(data.w.max.snp.genos,  age, by = "row.names")
  rownames(data.w.max.snp.genos) <- data.w.max.snp.genos[,"Row.names"]
  data.w.max.snp.genos <- data.w.max.snp.genos[,-1]
  
  colnames(data.w.max.snp.genos) <- c("expr", "genotype", "sex", "age")
  
  data.w.max.snp.genos$sex <- as.factor(data.w.max.snp.genos$sex)
  data.w.max.snp.genos$age <- as.factor(data.w.max.snp.genos$age)
  data.w.max.snp.genos$genotype <- as.factor(data.w.max.snp.genos$genotype)
  
  # Make expression vs age by genotype plots.
  par(mfrow = c(3,1))
  
  position_label <- paste0("by Genotype at Chr", snpinfo$chr , ", ", snpinfo$pos*10^6, "bp")
  # Both sexes
  summary_df <- summarySE(data.w.max.snp.genos, measurevar="expr", groupvars=c("genotype", "age"), na.rm = TRUE)
  summary_df <- na.omit(summary_df)
  p <- ggplot(summary_df, aes(x = age,  y = expr, color = genotype, group = genotype)) +
    geom_errorbar(aes(ymin=expr-se, ymax=expr+se), width=.1, position=position_dodge(0.1)) +
    geom_line(position=position_dodge(0.1)) +
    geom_point(position=position_dodge(0.1), size=3) +
    ggtitle(label = paste0(symbol, " (", id, ") ", "Expression vs Age \n", position_label)) +
    ylab(label = "Mean Normalized Expression Level") +
    xlab(label = "Age (months)") +
    labs(color = "Genotype") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  
  data.w.max.snp.genos.males <- data.w.max.snp.genos %>% filter(sex == "M")
  summary_df <- summarySE(data.w.max.snp.genos.males, measurevar="expr", groupvars=c("genotype", "age"), na.rm = TRUE)
  summary_df <- na.omit(summary_df)
  p <- ggplot(summary_df, aes(x = age,  y = expr, color = genotype, group = genotype)) +
    geom_errorbar(aes(ymin=expr-se, ymax=expr+se), width=.1, position=position_dodge(0.1)) +
    geom_line(position=position_dodge(0.1)) +
    geom_point(position=position_dodge(0.1), size=3) +
    ggtitle(label = paste0(symbol, " (", id, ") ", "Expression vs Age (Males) \n", position_label)) +
    ylab(label = "Mean Normalized Expression Level") +
    xlab(label = "Age (months)") +
    labs(color = "Genotype") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  
  data.w.max.snp.genos.females <- data.w.max.snp.genos %>% filter(sex == "F")
  summary_df <- summarySE(data.w.max.snp.genos.females, measurevar="expr", groupvars=c("genotype", "age"), na.rm = TRUE)
  summary_df <- na.omit(summary_df)
  p <- ggplot(summary_df, aes(x = age,  y = expr, color = genotype, group = genotype)) +
    geom_errorbar(aes(ymin=expr-se, ymax=expr+se), width=.1, position=position_dodge(0.1)) +
    geom_line(position=position_dodge(0.1)) +
    geom_point(position=position_dodge(0.1), size=3) +
    ggtitle(label = paste0(symbol, " (", id, ") ", "Expression vs Age (Females) \n", position_label)) +
    ylab(label = "Mean Normalized Expression Level") +
    xlab(label = "Age (months)") +
    labs(color = "Genotype") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
}
##############################################
##############################################
dir.create(output_dir)
setwd(output_dir)

hotspot_output_dirs <- list.dirs()
hotspot_output_dirs <- hotspot_output_dirs[grepl("hotspots_output", hotspot_output_dirs)]
hotspot_output_dirs <- hotspot_output_dirs[!grepl("allele", hotspot_output_dirs)]
# Iterate over each output directory. For each output directory, iterate over all
# hotspots in it.

hotspot_output_dirs <- hotspot_output_dirs[grepl("protein", hotspot_output_dirs)]
foreach (hotspot_output_dir = hotspot_output_dirs) %dopar% {
  # for (hotspot_output_dir in hotspot_output_dirs) {
  # hotspot_output_dir <- "./protein_age_hotspots_output"
  setwd(hotspot_output_dir)
  # Get hotspot gene lists in this current directory 
  # (each gene list corresponds to a hotspot)
  gene_list_files <- list.files()
  gene_list_files <- gene_list_files[grepl("_hotspot_genes.csv" , gene_list_files)]
  
  dir.create("./allele_effect_outputs/")
  
  # Iterate over each set of gene lists in current directory (current test type)
  foreach (file_suffix = gene_list_files) %dopar% {
    # for (file_suffix in gene_list_files) {
    if (grepl("sex", file_suffix)) {
      intcovar_name <- "Sex"
    } else if (grepl("age", file_suffix)) {
      intcovar_name <- "Age"
    } else {
      intcovar_name <- "None"
    }
    
    if (grepl("protein", file_suffix)) {
      expr_type <- "protein"; qtl_type <- "pQTL"
    } else {
      expr_type <- "mrna"; qtl_type <- "eQTL"
    }
    
    hotspot_annots <- read.csv(paste0("./", file_suffix))
    hotspot_annots <- hotspot_annots[order(-hotspot_annots$qtl_lod),]
    # interesting_genes <- hotspot_annots[1:5,] # Only make plots for genes with the largest 5 qtls.
    interesting_genes <- hotspot_annots
    
    # Custom symbol search
    # interesting_genes <- hotspot_annots[hotspot_annots$gene_symbol == "Ttn",]
    # if (nrow(interesting_genes) == 0) {next}
    #
    
    setwd("./allele_effect_outputs/")
    # for (i in 1:nrow(interesting_genes)) {
    foreach (i = 1:nrow(interesting_genes)) %dopar% {
      row <- interesting_genes[i,]
      id <- row$lodcolumn
      symbol <- row$gene_symbol
      qtl_chr <- row$qtl_chr
      qtl_pos <- row$qtl_pos
      
      if (expr_type == "protein") {
        expr.data <-  expr.protein[,as.character(id)]
        symbol <- toupper(symbol)
      } else {
        expr.data <- expr.mrna[,as.character(id)]
      }

      # id <- "ENSMUSG00000054342"; symbol <- "Kcnn4"; qtl_pos <- 65.506169
      # intcovar_name <- "Age"; expr_type <- "mrna"; qtl_type <- "eQTL"; expr.data <- expr.mrna[,as.character(id)]

      pdf(paste0(symbol, "_story_", id ,"_", tolower(intcovar_name), ".pdf"))
      ###################################################################
      ###################################################################
      
      #### QTL scans and plots
      par(mfrow=c(2, 1))
      if (intcovar_name != "None") {
        intcovar <- eval(parse(text = tolower(intcovar_name)))
        scan_output_full <- scan1(genoprobs = genoprobs, pheno = expr.data, kinship = kin_mat_list, addcovar = cbind(sex, age), intcovar = intcovar, cores = num_cores)
        scan_output_add <- scan1(genoprobs = genoprobs, pheno = expr.data, kinship = kin_mat_list, addcovar = cbind(sex, age), cores = num_cores)
        scan_output_effect <- scan_output_full - scan_output_add
        
        # Plot
        plot(x = scan_output_full, map = gmap, col = color[1])
        plot(x = scan_output_add, map = gmap, col = color[2], add = TRUE)
        plot(x = scan_output_effect, map = gmap, col = color[3], add = TRUE)
        
        # Plot title
        max_pos <- max(scan_output_full, gmap)$pos
        if (max_pos > median(gmap[[1]])) {pos <- "topleft"} else {pos <- "topright"}
        legend(pos, lwd = 2, col = color, c("Full", "Additive", "Interaction Effect"),  
               bg="gray90", lty=c(1,1,1), text.font = 10)
        title(paste0(symbol, " (", id, ") ", qtl_type, "\n Interactive Covariate = ", intcovar_name))
        
        qtl_chr <- as.numeric(max(scan_output_effect, gmap)$chr)
      } else { # intcovar_name == "None"
        
        scan_output_add <- scan1(genoprobs = genoprobs, 
                                 pheno = expr.data, 
                                 kinship = kin_mat_list, 
                                 addcovar = cbind(sex, age), 
                                 cores = num_cores)
        # Plot
        plot(x = scan_output_add, map = gmap, col = color[2])
        
        # Plot title
        max_pos <- max(scan_output_add, gmap)$pos
        if (max_pos > median(gmap[[1]])) {pos <- "topleft"} else {pos <- "topright"}
        title(paste0(symbol, " (", id, ") ", qtl_type))
        
        qtl_chr <- as.numeric(max(scan_output_add, gmap)$chr)
      }
      
      
      
      ##### Run allele effects. Run for all mice and then run stratified by our interactive covariate
      # (either age or sex or none. For the latter case we don't stratify)
      title <- paste0(symbol, " (", id, ") Expression")
      coef <- scan1blup(genoprobs = genoprobs[,qtl_chr],
                        pheno = expr.data,
                        kinship = kin_mat_list[[qtl_chr]],
                        addcovar = cbind(sex, age))
      plot_coefCC(x = coef, map = gmap[qtl_chr], main = title)
      
      # Stratify by our covariate.
      if (intcovar_name == "Sex") {    # Stratify on sex.
        
        # Get mouse ids for each sex group.
        annot.sample.M <- annot.sample %>% filter(annot.sample$Sex == "M")
        ids.M <- annot.sample.M$Mouse.ID
        ids.M <- ids.M[ids.M %in% rownames(kin_mat_list[[qtl_chr]])]
        
        annot.sample.F <- annot.sample %>% filter(annot.sample$Sex == "F")
        ids.F <- annot.sample.F$Mouse.ID
        ids.F <- ids.F[ids.F %in% rownames(kin_mat_list[[qtl_chr]])]
        
        par(mfrow = c(2,2))
        ########################
        # Make QTL and allele effect plots.
        ########################
        
        for (x in list(list(ids.M, "Males"), list(ids.F, "Females"))) {
          
          group_ids <- x[[1]]
          group_name <- x[[2]]
          
          ###################################
          # QTL scan
          scan_output_add <- scan1(genoprobs = genoprobs[,qtl_chr][group_ids], 
                                   pheno = expr.data[group_ids],
                                   kinship = kin_mat_list[[qtl_chr]][group_ids, group_ids], 
                                   addcovar = age[group_ids], cores = num_cores)
          plot(x = scan_output_add, map = gmap[qtl_chr], col = color[2])
          
          title(paste0(symbol, " (", id, ") ", qtl_type,"\n", group_name))
          
          ###################################
          # Allele effects
          title <- paste0(symbol, " (", id, "): ", group_name)
          coef <- scan1blup(genoprobs =  genoprobs[,qtl_chr][group_ids],
                            pheno = expr.data[group_ids],
                            kinship = kin_mat_list[[qtl_chr]][group_ids, group_ids],
                            addcovar = age[group_ids], cores = num_cores)
          plot_coefCC(x = coef, map = gmap[qtl_chr], main = title)
        } # for (x in list(c(ids.6, "6 Months"), c(ids.12, "12 Months"), c(ids.18, "18 Months"))) 
        
        
      } else if (intcovar_name == "Age") { # stratify by age
        
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
        
        par(mfrow = c(3,2))
        
        ########################
        # Make QTL and allele effect plots.
        ########################
        
        for (x in list(list(ids.6, "6 Months"), list(ids.12, "12 Months"), list(ids.18, "18 Months"))) {
          
          group_ids <- x[[1]]
          group_name <- x[[2]]
          
          ###################################
          # QTL scan
          scan_output_add <- scan1(genoprobs = genoprobs[,qtl_chr][group_ids], 
                                   pheno = expr.data[group_ids],
                                   kinship = kin_mat_list[[qtl_chr]][group_ids, group_ids], 
                                   addcovar = sex[group_ids], cores = num_cores)
          plot(x = scan_output_add, map = gmap[qtl_chr], col = color[2])
          
          max_pos <- max(scan_output_full, gmap[qtl_chr])$pos
          if (max_pos > median(gmap[qtl_chr][[1]])) {pos <- "topleft"} else {pos <- "topright"}
          title(paste0(symbol, " (", id, ") ", qtl_type,"\n", group_name))
          
          ###################################
          # Allele effects
          title <- paste0(symbol, " (", id, "): ", group_name)
          coef <- scan1blup(genoprobs =  genoprobs[,qtl_chr][group_ids],
                            pheno = expr.data[group_ids],
                            kinship = kin_mat_list[[qtl_chr]][group_ids, group_ids],
                            addcovar = sex[group_ids], cores = num_cores)
          plot_coefCC(x = coef, map = gmap[qtl_chr], main = title)
        } # for (x in list(c(ids.6, "6 Months"), c(ids.12, "12 Months"), c(ids.18, "18 Months"))) 
      } # if (intcovar_name == "age") else
      
      ###################################################################
      ###################################################################
      
      
      # Do association mapping (impute founder genotypes onto each DO mouse's genome)
      #########
      mapping_start = qtl_pos - 2
      mapping_end = qtl_pos + 2
      
      # Get genes under the peak
      query_genes = create_gene_query_func(dbfile = mouse_genes.sqlite, filter = "source='MGI'")
      genes = query_genes(qtl_chr, mapping_start, mapping_end)
      
      par(mfrow = c(1,1))
      if (intcovar_name != "None") {
        assoc_c <- scan1snps(genoprobs = genoprobs[,qtl_chr], map = gmap[qtl_chr],
                             pheno = expr.data, kinship =  kin_mat_list[[qtl_chr]],
                             addcovar = cbind(sex, age), intcovar = intcovar,
                             query_func = query_func, chr = qtl_chr,
                             start = mapping_start, end = mapping_end, keep_all_snps = TRUE, cores = num_cores)
        
        # Plot Manhattan plot + genes nearby
        plot_snpasso(scan1output = assoc_c$lod, snpinfo = assoc_c$snpinfo, genes = genes,
                     main = paste0(symbol, " (", id ,") Expression \n Interactive Covariate = ", intcovar_name))
        
      } else { # intcovar_name == "None")
        assoc_c <- scan1snps(genoprobs = genoprobs[,qtl_chr], map = gmap[qtl_chr],
                             pheno = expr.data, kinship =  kin_mat_list[[qtl_chr]],
                             addcovar = cbind(sex, age),
                             query_func = query_func, chr = qtl_chr,
                             start = mapping_start, end = mapping_end, keep_all_snps = TRUE, cores = num_cores)
        
        # Plot Manhattan plot + genes nearby
        plot_snpasso(scan1output = assoc_c$lod, snpinfo = assoc_c$snpinfo, genes = genes,
                     main = paste0(symbol,  " (", id ,") Expression"))
      }
      
      # Get SNP with highest LOD.
      max_snp = rownames(assoc_c$lod)[which.max(assoc_c$lod[,1])]
      
      # Get SNP info for that SNP.
      snpinfo = index_snps(gmap[qtl_chr], assoc_c$snpinfo[assoc_c$snpinfo$snp_id == max_snp,])
      
      # Get SNP probs. (Convert haplotype probabilities to allele probabilities (Main vs minor))
      snpprob = genoprob_to_snpprob(genoprobs[,qtl_chr], snpinfo = snpinfo)
      # This give you two columns, A & B. You should be able to use either one for a genotype by phenotype plot.
      
      # Round genotypes to 0, .5, or 1
      anchors <- c(0, .5, 1)
      for (k in 1:nrow(snpprob[[1]])) {
        genotype_prob_A <- snpprob[[1]][k,1,]
        genotype_prob_B <- snpprob[[1]][k,2,]
        snpprob[[1]][k,1,] <- anchors[which.min(abs(anchors - genotype_prob_A))]
        snpprob[[1]][k,2,] <- anchors[which.min(abs(anchors - genotype_prob_B))]
      } 
      
      # Convert to As and Bs
      genotypes_snp <- as.data.frame(snpprob[[1]])
      genotypes_snp[,"genotype"] <- NA
      for (k in 1:nrow(genotypes_snp)) {
        if (genotypes_snp[k,1] == 1) {
          genotypes_snp[k,"genotype"] <- "AA"
        } else if (genotypes_snp[k,1] == 0.5)  {
          genotypes_snp[k,"genotype"] <- "AB"
        } else {
          genotypes_snp[k,"genotype"] <- "BB"
        }
      }
      genotypes_snp <- genotypes_snp %>% select('genotype')
      genotypes_snp$genotype <- as.factor(genotypes_snp$genotype)
      
      # Make genotype plots
      expr_v_age_by_genotype_plots(expr.data = expr.data,  symbol = symbol, id = id, 
                                   snpinfo = snpinfo, genotypes_snp = genotypes_snp, annot.sample = annot.sample) 
      
      # If we are looking at protein expression, also
      # make mrna expression v age by genotype plots
      # and then do an eQTL scan
      # (If we are looking at mrna expression, we can't make a corresponding protein expr v age plot
      # since each mrna can map to several proteins)
      
      if (expr_type == "protein") {
        gene_id <- row$gene_id
        expr.data.mrna <- expr.mrna[,row$gene_id]
        
        expr_v_age_by_genotype_plots(expr.data = expr.data.mrna,  symbol = tools::toTitleCase(tolower(symbol)), id = gene_id, 
                                     snpinfo = snpinfo, genotypes_snp = genotypes_snp, annot.sample = annot.sample) 
        
        par(mfrow = c(2,1))
        if (intcovar_name != "None") {
          
          intcovar <- eval(parse(text = tolower(intcovar_name)))
          scan_output_full <- scan1(genoprobs = genoprobs[,qtl_chr], pheno = expr.data, kinship = kin_mat_list[[qtl_chr]], addcovar = cbind(sex, age), intcovar = intcovar, cores = num_cores)
          scan_output_add <- scan1(genoprobs = genoprobs[,qtl_chr], pheno = expr.data, kinship = kin_mat_list[[qtl_chr]], addcovar = cbind(sex, age), cores = num_cores)
          scan_output_effect <- scan_output_full - scan_output_add
          
          # Plot
          plot(x = scan_output_full, map = gmap[qtl_chr], col = color[1])
          plot(x = scan_output_add, map = gmap[qtl_chr], col = color[2], add = TRUE)
          plot(x = scan_output_effect, map = gmap[qtl_chr], col = color[3], add = TRUE)
          
          # Plot title
          max_pos <- max(scan_output_full, gmap[qtl_chr])$pos
          if (max_pos > median(gmap[qtl_chr][[1]])) {pos <- "topleft"} else {pos <- "topright"}
          legend(pos, lwd = 2, col = color, c("Full", "Additive", "Interaction Effect"))
          title(paste0(symbol, " (", id, ") pQTL \n Interactive Covariate = ", intcovar_name))
          
          #####################
          intcovar <- eval(parse(text = tolower(intcovar_name)))
          scan_output_full <- scan1(genoprobs = genoprobs[,qtl_chr], pheno = expr.data.mrna, kinship = kin_mat_list[[qtl_chr]], addcovar = cbind(sex, age), intcovar = intcovar, cores = num_cores)
          scan_output_add <- scan1(genoprobs = genoprobs[,qtl_chr], pheno = expr.data.mrna, kinship = kin_mat_list[[qtl_chr]], addcovar = cbind(sex, age), cores = num_cores)
          scan_output_effect <- scan_output_full - scan_output_add
          
          # Plot
          plot(x = scan_output_full, map = gmap[qtl_chr], col = color[1])
          plot(x = scan_output_add, map = gmap[qtl_chr], col = color[2], add = TRUE)
          plot(x = scan_output_effect, map = gmap[qtl_chr], col = color[3], add = TRUE)
          
          # Plot title
          max_pos <- max(scan_output_full, gmap[qtl_chr])$pos
          if (max_pos > median(gmap[qtl_chr][[1]])) {pos <- "topleft"} else {pos <- "topright"}
          legend(pos, lwd = 2, col = color, c("Full", "Additive", "Interaction Effect"))
          title(paste0(tools::toTitleCase(tolower(symbol)), " (", gene_id, ") eQTL \n Interactive Covariate = ", intcovar_name))
          
        } else { # intcovar_name == "None"
          scan_output_add <- scan1(genoprobs = genoprobs[,qtl_chr], 
                                   pheno = expr.data, 
                                   kinship = kin_mat_list[[qtl_chr]], 
                                   addcovar = cbind(sex, age), 
                                   cores = num_cores)
          # Plot
          plot(x = scan_output_add, map = gmap[qtl_chr], col = color[2])
          
          # Plot title
          title(paste0(symbol, " (", id, ") pQTL \n (Additive)"))
          
          #######
          scan_output_add <- scan1(genoprobs = genoprobs[,qtl_chr], 
                                   pheno = expr.data.mrna, 
                                   kinship = kin_mat_list[[qtl_chr]], 
                                   addcovar = cbind(sex, age), 
                                   cores = num_cores)
          # Plot
          plot(x = scan_output_add, map = gmap[qtl_chr], col = color[2])
          
          # Plot title
          title(paste0(tools::toTitleCase(tolower(symbol)), " (", gene_id, ") eQTL \n (Additive)"))
        }
      } #  if (expr_type == "protein") {
      
      dev.off()
    } # for (i in 1:nrow(interesting_genes)) {
    setwd("..")
  } # for (file_suffix in files)
  setwd("..") # Go back up two directories
} # for (hotspot_output_dir in hotspot_output_dirs){
