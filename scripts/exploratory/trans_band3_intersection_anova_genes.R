low_low <- read.csv("~/do_heart/results/anova/age_corrs_mrna_low_protein_low_genes.csv")
high_high<- read.csv("~/do_heart/results/anova/age_corrs_mrna_high_protein_high_genes.csv")
high_low <- read.csv("~/do_heart/results/anova/age_corrs_mrna_high_protein_low_genes.csv")
low_high <- read.csv("~/do_heart/results/anova/age_corrs_mrna_low_protein_high_genes.csv")

protein_age_trans_chr3 <- read.csv("~/do_heart/results/hotspots/protein_age_hotspots_output/protein_age_trans_chr3_hotspot_genes.csv")
intersection_low_low <- protein_age_trans_chr3[protein_age_trans_chr3$gene_id %in% low_low$gene_id,]
intersection_low_low$gene_symbol
# [1] Myh6   Ttn    Myh6   Myh7b  Obscn  Ank2   Ttn    Map7d1 Ttn    Ttn    Obscn 
intersection_low_low$gene_id

intersection_high_high <- protein_age_trans_chr3[protein_age_trans_chr3$gene_id %in% high_high$gene_id,]
intersection_high_high$gene_symbol
# Nothing

intersection_high_low <- protein_age_trans_chr3[protein_age_trans_chr3$gene_id %in% high_low$gene_id,]
intersection_high_low$gene_symbol
# Ahsg

intersection_low_high <- protein_age_trans_chr3[protein_age_trans_chr3$gene_id %in% low_high$gene_id,]
intersection_low_high$gene_symbol
# Nothing
