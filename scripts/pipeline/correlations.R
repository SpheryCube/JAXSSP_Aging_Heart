# correlations
# Author: Daniel Alfonsetti
# Date: June 18, 2018

# Computes raw, sex adjusted, and age adjusted correlations between protein and mRNA expression levels.
# and plots the results

#######################
# Reset environment
rm(list=ls())

# Using local directories
# input <- "~/do_heart/data/cleaned_data.RData"
# output_prefix <- "~/do_heart/scripts/pipeline/output/" 

# Using HPC directories
input <- "/fastscratch/c-alfond/do_heart/data/cleaned_data.RData"
load(input)
output_prefix <- "/fastscratch/c-alfond/do_heart/results/"

#######################
library(ggplot2)
library(dplyr)
#######################

# Can only compute correlations for protein and mRNA that come in pairs.

# Intialize a dataframe to store our results.
corr.df <- data.frame(matrix(ncol=9, nrow = 0))
names(corr.df) <- c("gene.id", "protein.id", "symbol", "chr", "start", "end", "r.raw", "r.sex_adj", "r.age_adj")


# Iterate over genes row by row (only for mRNA and protein data that matches)
for (i in 1:N$complete) {
  
  gene_row <- annot.protein[i,]
  
  # Pull out relevant data
  current_symbol = gene_row$symbol
  current_chr = gene_row$chr
  current_gene_id <- gene_row$gene_id
  current_protein_id <- gene_row$protein_id # Get the correspoinding protein_id
  current_start <- gene_row$start
  current_end <- gene_row$end
  
  mrna_expr_data <- expr.mrna[, current_gene_id]
  protein_expr_data <- expr.protein[, current_protein_id]
  
  # Create gene-specific data frame holding data to send to our linear models.
  df <- cbind(annot.sample$Batch, annot.sample$Generation, annot.sample$Age, annot.sample$Sex, annot.sample$TMT, mrna_expr_data, protein_expr_data)
  df <- as.data.frame(df)
  
  colnames(df) <- c("Batch", "Generation", "Age", "Sex", "TMT", "mrna.expr", "protein.expr")
  
  # Get raw
  r.raw <- cor(mrna_expr_data, protein_expr_data)
  
  # Get sex adjusted correlation
  fit.prot.sex <- lm(protein_expr_data ~ Sex + Batch + Generation + TMT, data = df)
  res.prot.sex <- residuals(fit.prot.sex)
  
  fit.mrna.sex <- lm(mrna_expr_data ~ Sex + Batch + Generation + TMT, data = df)
  res.mrna.sex <- residuals(fit.mrna.sex)
  
  r.sex_adj <- cor(res.mrna.sex, res.prot.sex)
  
  # Get age adjusted correlation
  fit.prot.age <- lm(protein_expr_data ~ Age + Batch + Generation + TMT, data = df)
  res.prot.age <- residuals(fit.prot.age)
  
  fit.mrna.age <- lm(mrna_expr_data ~ Age + Batch + Generation + TMT, data = df)
  res.mrna.age <- residuals(fit.mrna.age)
  
  r.age_adj <- cor(res.mrna.age, res.prot.age)
  
  
  # Store results in the 'results' dataframe we are building up row by row
  corr.df[nrow(corr.df)+1,] <- list(current_gene_id, current_protein_id, current_symbol, current_chr, current_start, current_end, r.raw, r.sex_adj, r.age_adj)
  
}

# Save the results dataframe
write.csv(corr.df, file=paste0(output_prefix, "correlations_output.csv"), row.names = FALSE)

# Diagnostics
# for (colname in c("r.raw", "r.sex_adj", "r.age_adj"))
# {
#   title <- paste(colname, "correlation histogram (for complete data)", sep = " ")
#   plot <- ggplot(data= corr.df, aes(x = eval(parse(text= colname)))) +
#     geom_histogram(bins = 100,
#                    col="blue",
#                    fill="green",
#                    alpha = .2) +
#     ggtitle(title) +
#     xlab(colname)
#   print(plot)
# }
