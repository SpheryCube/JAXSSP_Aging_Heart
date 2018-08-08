# Author: Daniel Alfonsetti
# Date: July 12, 2018 
# Description: Follow up to the 'mRNA_protein_cors.R' script.
# Does the overall correlation between mrna and their respective transcripts change with age? 
# Conclusion: Yes, but its sex dependent.

# Reference(s): https://stats.stackexchange.com/questions/171496/r-test-for-correlation-with-a-covariate
#######################
# Reset environment
rm(list=ls())

# Using local directories
input <- "~/do_heart/data/cleaned_data.RData"
load(input)
output_dir <- "~/do_heart/results/mrna_protein_cors_by_age_sex/"
setwd(output_dir)
# Using HPC directories
# input <- "/home/c-alfond/do_heart/data/cleaned_data.RData"
# load(input)
# output_prefix <- "/home/c-alfond/do_heart/scripts/pipeline/output/"

#######################
library(ggplot2)
library(dplyr)
library(broom)
library(grid)
library(ppcor)
library(RColorBrewer)
#######################
#######################

# Compute the correlation within age groups, draw box plots (with jittered points overlay),
# fit the regression on age and evaluate significance of the slope.

# Intialize a dataframe to store our results.
corr.df <- data.frame(matrix(ncol=24, nrow = 0))
names(corr.df) <- c("gene.id", "protein.id", "symbol", 
                    "chr", "start", "end", 
                    "r.M.6", "r.M.12", "r.M.18", "d6.12.M", "d6.18.M", "d12_18.M", 
                    "r.F.6", "r.F.12", "r.F.18","d6.12.F", "d6.18.F", "d12.18.F",
                    "r.6", "r.12","r.18", "d6.12", "d6.18", "d12.18")

# Can only compute correlations for protein and mRNA that come in pairs.
# Iterate over genes row by row 
for (i in 1:N$complete) { #1:N$pairs. Can't use if using pcor function.
  gene_row <- annot.protein[i,]
  
  # Pull out relevant data
  current_symbol = gene_row$symbol
  current_chr = as.character(gene_row$chr)
  current_gene_id <- gene_row$gene_id
  current_protein_id <- gene_row$protein_id # Get the correspoinding protein_id
  current_start <- gene_row$start
  current_end <- gene_row$end
  

  #######################
  #######################
  # Get protein and transcript expression data for current gene. 
  # bind annotation data to it so we have covariates to run partial correlations with.
  exprM <- expr.mrna[, current_gene_id]
  exprP <- expr.protein[, current_protein_id]
  mrna_expr_data <- cbind(exprM, annot.sample)
  protein_expr_data <- cbind(exprP, annot.sample)
  expr_data_w_annots <- cbind(exprM, exprP, 
               as.numeric(as.factor(annot.sample[,'TMT'])), 
               as.numeric(as.factor(annot.sample[,'Batch'])), 
               as.numeric(as.factor(annot.sample[,'Generation'])))
  expr_data_w_annots <- as.matrix(expr_data_w_annots)
  #######################
  # Separate groups by age and sex
  for (sex in c("M", "F")) {
    for (age in c(6, 12, 18)) {

      filtered_mice <- annot.sample %>% filter(Sex == sex, Age == age)
      
      expr_data_w_annots_filtered <- expr_data_w_annots[filtered_mice$Mouse.ID,]
      
      pcor_out <- pcor(expr_data_w_annots_filtered)

      result <- pcor_out$estimate[1,2]
      assign(paste("r", sex, as.character(age), sep = "."), result)
    }
  }
  #######################
  # Now only seperate on age, not sex.
  expr_data_w_annots <- cbind(expr_data_w_annots, as.numeric(as.factor(annot.sample[,'Sex'])))
  for (age in c(6, 12, 18)) {
    filtered_mice <- annot.sample %>% filter(Age == age)
    expr_data_w_annots_filtered <- expr_data_w_annots[filtered_mice$Mouse.ID,]
    pcor_out <- pcor(expr_data_w_annots_filtered) # Calculate partial correlation coefficients 
    result <- pcor_out$estimate[1,2]
    assign(paste("r", as.character(age), sep = "."), result)
  }
  
  #######################
  #######################
  # Calculate correlation differences
  d6.12.M <- r.M.12 - r.M.6; 
  d6.18.M <- r.M.18 - r.M.6; 
  d12.18.M <- r.M.18 - r.M.12
  
  d6.12.F <- r.F.12 - r.F.6; 
  d6.18.F <- r.F.18 - r.F.6; 
  d12.18.F <- r.F.18 - r.F.12
  
  d6.12 <- r.12 - r.6; 
  d6.18 <- r.18 - r.6; 
  d12.18 <- r.18 - r.12
  
  ###########
  # Store results in the 'results' dataframe we are building up row by row
  corr.df[nrow(corr.df)+1,] <- list(current_gene_id, current_protein_id, current_symbol, 
                                    current_chr, current_start, current_end, 
                                    r.M.6, r.M.12, r.M.18, d6.12.M, d6.18.M, d12.18.M,
                                    r.F.6, r.F.12, r.F.18, d6.12.F, d6.18.F, d12.18.F,
                                    r.6, r.12, r.18, d6.12, d6.18, d12.18)
  
}

# Save the results dataframe
write.csv(corr.df, file=paste0(output_dir, "mrna_protein_cors_by_age_sex.csv"), row.names = FALSE)

########################################
# Graphing
Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                 cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*", " "))



########################################
par(mfrow=c(3,1))

pdf("mRNA_protein_cors_by_age_boxplots.pdf")

# Male correlations by age
corr.df.male <- corr.df[,c("r.M.6", "r.M.12", "r.M.18")]
corr.df.male <- stack(corr.df.male)

age_group <- as.numeric(corr.df.male$ind) # Levels: 1, 2, 3
lm_male <- lm(values ~ age_group, data = corr.df.male)
result_male <- tidy(lm_male)
result_male
# term     estimate   std.error statistic       p.value
# 1 (Intercept)  0.233824331 0.008007197 29.201770 1.716012e-177
# 2   age_group -0.008932868 0.003706612 -2.409982  1.597763e-02

Signif <- symnum(result_male$p.value[2], corr = FALSE, na = FALSE,
                 cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*", " "))
annot_sig <- grobTree(textGrob(Signif, x= unit(0.90, "npc"), y = unit(0.96, "npc"),
                               gp=gpar(col='red', fontsize=30, fontface="italic")))
annot_m <- grobTree(textGrob(paste0("Age coefficient: ", round(result_male$estimate[2], 3)), x= unit(0.70, "npc"), y = unit(0.97, "npc"),
                             gp=gpar(col='red', fontsize=15, fontface="italic")))

ggplot(corr.df.male, aes(x = ind, y = values)) +
  scale_color_brewer(palette="Dark2") +
  geom_jitter(alpha = 0.5, aes(color = ind)) +
  geom_boxplot(outlier.shape = NA, color = "navyblue", alpha = 0, width = 0.4) + 
  geom_violin(alpha = 0) +
  scale_x_discrete(name = "Age Groups", labels = c("6 Months", "12 Months", "18 Months")) + 
  ylab("mRNA-Protein Correlations") + 
  theme_classic() +
  theme(legend.position="none") +
  ggtitle("mRNA-Protein correlations for males grouped by age") +
  annotation_custom(annot_sig) +
  annotation_custom(annot_m)
  
print(paste0("Age slope significant for males: ", result_male$p.value[2]))


##########
##########
# Female correlations by age

corr.df.female <- corr.df[c("r.F.6", "r.F.12", "r.F.18")]
corr.df.female <- stack(corr.df.female)

age_group <- as.numeric(corr.df.female$ind) # Levels: 1, 2, 3
lm_female <- lm(values~ age_group, data = corr.df.female)
result_female <- tidy(lm_female)
result_female
# term    estimate   std.error statistic       p.value
# 1 (Intercept)  0.27719071 0.008152395  34.00114 1.575111e-235
# 2   age_group -0.02993927 0.003773825  -7.93340  2.447158e-15

Signif <- symnum(result_female$p.value[2], corr = FALSE, na = FALSE,
                 cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*", " "))
annot_sig <- grobTree(textGrob(Signif, x= unit(0.90, "npc"), y = unit(0.96, "npc"),
                               gp=gpar(col='red', fontsize=30, fontface="italic")))
annot_m <- grobTree(textGrob(paste0("Age coefficient: ", round(result_female$estimate[2], 3)), x= unit(0.70, "npc"), y = unit(0.97, "npc"),
                             gp=gpar(col='red', fontsize=15, fontface="italic")))

ggplot(corr.df.female, aes(x = ind, y = values)) +
  scale_color_brewer(palette="Dark2") +
  geom_jitter(alpha = 0.5, aes(color = ind)) +
  geom_boxplot(outlier.shape = NA, color = "navyblue", alpha = 0, width = 0.4) + 
  geom_violin(alpha = 0) +
  scale_x_discrete(name = "Age Groups", labels = c("6 Months", "12 Months", "18 Months")) + 
  ylab("mRNA-Protein Correlations") + 
  theme_classic() +
  theme(legend.position="none") +
  ggtitle("mRNA-Protein correlations for females grouped by age") +
  annotation_custom(annot_sig) +
  annotation_custom(annot_m)
##########
##########
# All correlations by age.

corr.df.age <- corr.df[c("r.6", "r.12", "r.18")]
corr.df.age <- stack(corr.df.age)

age_group <- as.numeric(corr.df.age$ind) # Levels: 1, 2, 3
lm_age <- lm(values~ age_group, data = corr.df.age)
result_age <- tidy(lm_age)
result_age
# term   estimate   std.error statistic       p.value
# 1 (Intercept)  0.2567527 0.007050322 36.417164 3.537868e-267
# 2   age_group -0.0188097 0.003263665 -5.763368  8.577448e-09

Signif <- symnum(result_age$p.value[2], corr = FALSE, na = FALSE,
                 cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                 symbols = c("***", "**", "*", " "))
annot_sig <- grobTree(textGrob(Signif, x= unit(0.90, "npc"), y = unit(0.96, "npc"),
                               gp=gpar(col='red', fontsize=30, fontface="italic")))
annot_m <- grobTree(textGrob(paste0("Age coefficient: ", round(result_age$estimate[2], 3)), x= unit(0.70, "npc"), y = unit(0.97, "npc"),
                             gp=gpar(col='red', fontsize=15, fontface="italic")))

ggplot(corr.df.age, aes(x = ind, y = values)) +
  scale_color_brewer(palette="Dark2") +
  geom_jitter(alpha = 0.5, aes(color = ind)) +
  geom_boxplot(outlier.shape = NA, color = "navyblue", alpha = 0, width = 0.4) + 
  geom_violin(alpha = 0) +
  scale_x_discrete(name = "Age Groups", labels = c("6 Months", "12 Months", "18 Months")) + 
  ylab("mRNA-Protein Correlations") + 
  theme_classic() +
  theme(legend.position="none") +
  ggtitle("mRNA-Protein correlations grouped by age") +
  annotation_custom(annot_sig) +
  annotation_custom(annot_m)


dev.off()


# Conclusion: Each summary shows a significant non-zero negative slope indicating that correlations between mRNA transcripts and 
# their respective proteins decrease in the heart as mice age.
##############
##############
# Make histograms of correlations by age and sex.

corr.df <- read.csv("mRNA_protein_cors_by_age.csv")

# Rearrange data to make histogram facet.
corrs_only <- corr.df[, c("r.M.6", "r.M.12", "r.M.18", "r.F.6", "r.F.12", "r.F.18")]
corrs_by_group <- stack(corrs_only)
corrs_by_group$ind <- as.character(corrs_by_group$ind)
for (i in 1:nrow(corrs_by_group)) {
  ind <- corrs_by_group[i,]$ind
  if (ind == "r.M.6") {corrs_by_group[i, "Age"] <- "6 Months"; corrs_by_group[i, "Sex"] <- "Male"}
  else if (ind == "r.M.12") {corrs_by_group[i, "Age"] <- "12 Months"; corrs_by_group[i, "Sex"] <- "Male"}
  else if (ind == "r.M.18") {corrs_by_group[i, "Age"] <- "18 Months"; corrs_by_group[i, "Sex"] <- "Male"}
  else if (ind == "r.F.6") {corrs_by_group[i, "Age"] <- "6 Months"; corrs_by_group[i, "Sex"] <- "Female"}
  else if (ind == "r.F.12") {corrs_by_group[i, "Age"] <- "12 Months"; corrs_by_group[i, "Sex"] <- "Female"}
  else {corrs_by_group[i, "Age"] <- "18 Months"; corrs_by_group[i, "Sex"] <- "Female"}
}
corrs_by_group$Age <- factor(corrs_by_group$Age, levels = c("6 Months", "12 Months", "18 Months"))

pdf("mrna_protein_cors_age_sex.pdf")
p <- ggplot(data = corrs_by_group, aes(values)) +
  geom_histogram(breaks=seq(-.75, 1, by=0.01),
                 aes(fill=..count..)) +
  theme_linedraw() + 
  facet_grid(vars(Age), vars(Sex)) +
  ggtitle("Protein-Transcript Correlations Grouped by Age and Sex") +
  xlab("Correlation Coefficients") +
  ylab("Count")
p
dev.off()

#################################
#################################
# Now make something like
#     mean_corr.M     sd_corr.M       se_corr.M     mean_corr.F  sd_corr.F    se_corr.F
# 6     ...              ... 
# 12    ...      
# 18    ...


corr_mean_sd_se.df <- data.frame(matrix(nrow = 3, ncol = 6))
colnames(corr_mean_sd_se.df) <- c("mean_corr.M", "sd_corr.M", "se_corr.M",
                                  "mean_corr.F", "sd_corr.F", "se_corr.F")
rownames(corr_mean_sd_se.df) <- c("6 Months", "12 Months", "18 Months")


for (age in c("6 Months", "12 Months", "18 Months")) {
   group_male <- corrs_by_group %>% filter(Age == age, Sex == "Male")
   M_mean <- mean(group_male$values)
   M_sd <- sd(group_male$values)
   M_sem <- sd(group_male$values)/sqrt(length(group_male$values))
   group_female <- corrs_by_group %>% filter(Age == age, Sex == "Female")
   F_mean <- mean(group_female$values)
   F_sd <- sd(group_female$values)
   F_sem <- sd(group_female$values)/sqrt(length(group_female$values))
   x <- c(M_mean, M_sd, M_sem, F_mean, F_sd, F_sem)
   x <- lapply(x, round, digits=5)
   corr_mean_sd_se.df[age,] <- x
}
write.csv(corr_mean_sd_se.df, "mrna_protein_cors_MeanSdSe_by_age_sex.csv")
corr_mean_sd_se.df
#             mean_corr.M sd_corr.M se_corr.M mean_corr.F sd_corr.F se_corr.F
# 6 Months      0.20985   0.26990   0.00544     0.25352   0.26682   0.00538
# 12 Months     0.24604   0.26435   0.00533     0.20478   0.26962   0.00544
# 18 Months     0.19198   0.24217   0.00488     0.19364   0.25695   0.00518

