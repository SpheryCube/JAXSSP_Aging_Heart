# Author: Daniel Alfonsetti (daniel.alfonsetti@gmail.com)
# Date: June 21, 2018
# Cleans data for analysis.
# Extended on from Gary Churchill's script. Now cleans genoprob data to ready for QTL analysis.

#######################
# # # Reset environment
# rm(list=ls())
# 
# load("~/do_heart/data/DO189_heart_v2_noprobs.RData") # # Load expression data
# load("~/do_heart/data/JAC_crosssectional_genoprobs_20180618.Rdata") # # Load genotype data
# 
# output_prefix <- "~/do_heart/data/" # # Place to put cleaned data.

############## (Reset environment using HPC)
rm(list=ls())

# Load expression data
load("/home/c-alfond/do_heart/data/DO189_heart_v2_noprobs.RData")
# Load genotype data
load("/home/c-alfond/do_heart/data/genotypes_crosssectional/JAC_crosssectional_genoprobs_20180618.Rdata")

# Place to put cleaned data.
output_prefix <- "/home/c-alfond/do_heart/data/"

#######################
library(dplyr)
library(qtl2)
#######################

# Convert chr to an ordered factor
annot.mrna$chr <- factor(annot.mrna$chr,
                         levels=c(as.character(1:19),"X","Y","MT"))
annot.protein$chr <- factor(annot.protein$ch,
                            levels=c(as.character(1:19),"X","Y","MT"))

# Change some names 
annot.protein <- dplyr::rename(annot.protein, protein_id = id)
annot.protein <- dplyr::rename(annot.protein, middle=middle_point)
annot.mrna <- dplyr::rename(annot.mrna, gene_id = id)
annot.mrna <- dplyr::rename(annot.mrna, middle=middle_point)
annot.samples <- dplyr::rename(annot.samples, TMT=Tag, Sample.ID=Sample.Number)

annot.sample <- transform(annot.samples, Batch=factor(substring(Batch,1,4)),
                          NYGC.ID=substring(NYGC.ID,first=11),
                          Sample.ID = as.integer(Sample.ID),
                          Age=as.integer(Age), 
                          Sex=factor(Sex),
                          Generation=factor(Generation,levels=c("G8","G9","G10","G11","G12")))


# # # Save sex and age numeric values
# Sex <- as.numeric(as.factor(annot.sample$Sex))
# names(Sex) <- rownames(annot.sample)
# Age <- as.numeric(annot.sample$Age)
# names(Age) <- rownames(annot.sample)

rm(annot.samples)


# Convert - to . in expr.mrna row names (mouse ids) so that IDs in genoprobs (probs) and expr.mrna match
names <- rownames(expr.mrna)
new_names <- gsub("-", ".", names)
rownames(expr.mrna) <- new_names

rownames(annot.sample) <- new_names

# Do the same for protein
names <- rownames(expr.protein)
new_names <- gsub("-", ".", names)
rownames(expr.protein) <- new_names

# Do the same for annnot.sample mouse ids
annot.sample$Mouse.ID <- new_names

# Rename map as gmap and probs as genoprobs
gmap <- map
genoprobs <- probs
rm(map)
rm(probs)

# For some reason, two mice have duplicated values.
# For each mouse, remove it's duplicate.
# Mice IDs "DO.0914" and "DO.0940" are repeated.
# In a given chromosome, data[177,,] == data[178,,]
# and data[202,,] == data[203,,].
for (chr in names(genoprobs)) {
  print(chr)
  data <- genoprobs[[chr]] # Get probabilities for this chromosome
  num_obs <- length(data[,1,1]) # Get how many observations there are.
  result <- unique(data[1:num_obs,,])  # remove duplicate observations
  genoprobs[[chr]] <- result # Save the result back into our genoprobs data structure
}
rm(data)
rm(num_obs)
rm(result)
rm(chr)

# Get kinship matrix
kin_mat_list <- calc_kinship(genoprobs, type = "loco", cores = 2)

# Save the cleaned genotype and expression data into new file.
save.image(paste0(output_prefix, "cleaned_data_with_genoprobs.RData"))

# Save a seperate version without genoprobs (genoprobs is a very large file and takes a long time to load)
rm(genoprobs)
rm(kin_mat_list)
save.image(paste0(output_prefix, "cleaned_data.RData"))


