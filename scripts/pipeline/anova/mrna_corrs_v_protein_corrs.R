# Author: Daniel Alfonsetti (daniel.alfonsetti@gmail.com)
# Date: July 26, 2018
# Description: Using the anova.R script output, plot mrna-age given sex correlations against
# protien-age given sex correlations. Color code genes that have significance correlations
# on both axes. For each group of genes (seperated by quadrant) that have significant correlations
# for both mrna-age and protein-age, perform enrichment analysis.
#######################
# Reset environment
rm(list = ls())
gc()
dev.off()
library(qtl2)
library(intermediate)
library(broom)
library(tidyverse)
library(grid)
library(allez)
library('org.Mm.eg.db')
library(clusterProfiler)
#######################

# Using local directories
input <- "~/do_heart/data/cleaned_data.RData"
load(input)
output_dir <- "~/do_heart/results/anova/"
dir.create(output_dir)
setwd(output_dir)

anovas <- read.csv("~/do_heart/results/anova_output.csv")
#######################
sex <- as.numeric(annot.sample$Sex) # 1 is female, 2 is male
names(sex) <- rownames(annot.sample)

age <- as.numeric(annot.sample$Age)
names(age) <- rownames(annot.sample)

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

df <- anovas %>% dplyr::select(protein_id, gene_id, symbol,  r.mRNA_Age.Sex, r.Prot_Age.Sex, p.mRNA_Age.Sex, p.Prot_Age.Sex)

df <- df %>% dplyr::mutate(sig = (df$p.mRNA_Age.Sex < 0.05 & df$p.Prot_Age.Sex < 0.05))

pdf("mrna_age_v_protein_age_corrs_given_sex.pdf")
plot <- ggplot(df, aes(x=r.mRNA_Age.Sex, y=r.Prot_Age.Sex, color = sig)) +
  geom_point(alpha = 0.3) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("mRNA-Age") +
  ylab("Protein-Age") +
  ggtitle("Protein-Age vs mRNA-Age Correlations given Sex") +
  xlim(-0.7, 0.7) +
  ylim(-1.1, 1.1)
print(plot)
dev.off()

#######################

sig_only <- df %>% filter(sig == TRUE)
mrna_low_protein_low <- sig_only %>% filter(r.Prot_Age.Sex < 0, r.mRNA_Age.Sex < 0)
mrna_low_protein_high <- sig_only %>% filter(r.Prot_Age.Sex > 0, r.mRNA_Age.Sex < 0)
mrna_high_protein_low <- sig_only %>% filter(r.Prot_Age.Sex < 0, r.mRNA_Age.Sex > 0)
mrna_high_protein_high <- sig_only %>% filter(r.Prot_Age.Sex > 0, r.mRNA_Age.Sex > 0)

write.csv(file = paste0(output_dir, "age_corrs_mrna_low_protein_low_genes.csv"), mrna_low_protein_low)
write.csv(file = paste0(output_dir, "age_corrs_mrna_low_protein_high_genes.csv"), mrna_low_protein_high)
write.csv(file = paste0(output_dir, "age_corrs_mrna_high_protein_low_genes.csv"), mrna_high_protein_low)
write.csv(file = paste0(output_dir, "age_corrs_mrna_high_protein_high_genes.csv"), mrna_high_protein_high)

# Do enrichment on each of these groups.
for (group_name in c("mrna_low_protein_low", "mrna_low_protein_high", "mrna_high_protein_low", "mrna_high_protein_high")) {
  
  group <- eval(parse(text = group_name))
    
  ###################################################
  # Cluster Profiler setup
  ###################################################
  
  # Create df that maps ENSEMBL, ENTREZ, and symbols to each other.
  # Contains all genes, not just significant genes.
  universe_df <- bitr(as.character(anovas$symbol),
                      fromType = "SYMBOL",
                      toType = c("ENTREZID", "ENSEMBL"),
                      OrgDb = org.Mm.eg.db)
  
  
  # Create df that contains only significant genes.
  sig_gene_df <- bitr(as.character(group$symbol),
                      fromType = "SYMBOL",
                      toType = c("ENTREZID", "ENSEMBL"),
                      OrgDb = org.Mm.eg.db)
  
  # ###########################################
  # clusterProfiler Part 1 - GO
  # ###########################################
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
  if (nrow(go_enrich_output) != 0)
  {
    write.csv(go_enrich_output, file = paste0(output_dir, group_name, "_age_gsea_go_clusterProfiler.csv"), row.names = FALSE)
  }

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
    write.csv(kegg_result, file = paste0(output_dir, group_name, "_age_gsea_kegg_clusterProfiler.csv"), row.names = FALSE)
  }
  # ##########################################
  # End of Allez enrichment analysis
  # ##########################################
  
} # (group in c(mrna_low_protein_low, mrna_low_protein_high, mrna_high_protein_low, mrna_high_protein_high)) {
