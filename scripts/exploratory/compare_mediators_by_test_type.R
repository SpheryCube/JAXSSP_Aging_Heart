# Daniel Alfonsetti
# August 7, 2018 
rm(list = ls())

library(tidyverse)
library(reshape2)
mrna_sex_mrna_mediators <- read.csv("~/do_heart/results/hotspots/mrna_sex_hotspots_output/mrna_sex_hotspot_mrna_mediators_r10.csv")
mrna_sex_protein_mediators <- read.csv("~/do_heart/results/hotspots/mrna_sex_hotspots_output/mrna_sex_hotspot_protein_mediators_r10.csv")

mrna_age_mrna_mediators <- read.csv("~/do_heart/results/hotspots/mrna_age_hotspots_output/mrna_age_hotspot_mrna_mediators_r10.csv")
mrna_age_protein_mediators <- read.csv("~/do_heart/results/hotspots/mrna_age_hotspots_output/mrna_age_hotspot_protein_mediators_r10.csv")


protein_sex_mrna_mediators <- read.csv("~/do_heart/results/hotspots/protein_sex_hotspots_output/protein_sex_hotspot_mrna_mediators_r10.csv")
protein_sex_protein_mediators <- read.csv("~/do_heart/results/hotspots/protein_sex_hotspots_output/protein_sex_hotspot_protein_mediators_r10.csv")

protein_age_mrna_mediators <- read.csv("~/do_heart/results/hotspots/protein_age_hotspots_output/protein_age_hotspot_mrna_mediators_r10.csv")
protein_age_protein_mediators <- read.csv("~/do_heart/results/hotspots/protein_age_hotspots_output/protein_age_hotspot_protein_mediators_r10.csv")


# Raw
l <- list(mrna_sex_mrna_mediators$diff, mrna_sex_protein_mediators$diff,
          mrna_age_mrna_mediators$diff, mrna_age_protein_mediators$diff,
          protein_sex_mrna_mediators$diff, protein_sex_protein_mediators$diff,
          protein_age_mrna_mediators$diff, protein_age_protein_mediators$diff)

# Ratioed on mean of interaction effect without mediator.
# l <- list(mrna_sex_mrna_mediators$diff/mean(mrna_sex_mrna_mediators$int_effect_lod), 
#           mrna_sex_protein_mediators$diff/mean(mrna_sex_protein_mediators$int_effect_lod),
#           mrna_age_mrna_mediators$diff/mean(mrna_age_mrna_mediators$int_effect_lod), 
#           mrna_age_protein_mediators$diff/mean(mrna_age_protein_mediators$int_effect_lod),
#           protein_sex_mrna_mediators$diff/mean(protein_sex_mrna_mediators$int_effect_lod), 
#           protein_sex_protein_mediators$diff/mean(protein_sex_protein_mediators$int_effect_lod),
#           protein_age_mrna_mediators$diff/mean(protein_age_mrna_mediators$int_effect_lod), 
#           protein_age_protein_mediators$diff/mean(protein_age_protein_mediators$int_effect_lod))



# https://stackoverflow.com/questions/38540949/r-repeatedly-cbind-a-matrix-a-vector-of-unequal-length-vector-goes-into-n
result <- lapply(l, function(x) x[1: max(sapply(l, length))]) %>% do.call(cbind, .) 

result <- as.data.frame(result)
colnames(result) <- c("mRNA_Sex: mRNA", "mRNA_Sex: Protein",
                      "mRNA_Age: mRNA", "mRNA_Age: Protein",
                      "Protein_Sex: mRNA", "Protein_Sex: Protein",
                      "Protein_Age: mRNA", "Protein_Age: Protein")
result <-  melt(result)

ggplot(data = result, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  theme_linedraw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Change in Interaction Effect with Mediator") +
  labs(fill = "OutcomeType_Intcovar: \n MediatorType") +
  ggtitle("Interactive QTL Mediation Analysis Result Summaries for Various Scan Types")

