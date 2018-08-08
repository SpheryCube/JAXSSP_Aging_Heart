##################################################################################
# R/3.4.4
# Originally created by Yuka Takemon; modified by Daniel Alfonsetti (July 13, 2018)

# The purpose of this script is to group categories from the allez GSEA output that were similiar.

# Load each file and create a matrix of Jaccard distances:
# J(A,B) = |A n B| / |A u B|
#        = count will give % similarity
# where A is the set of genes in one category and B is the set of genes in another category.

##################################################################################

# Reset environment
rm(list = ls())
set.seed(123)

# load libraries
library(tidyverse)
library(reshape2)
library(gplots)
library(foreach)

###
# set up directories (Local)
output_dir <- "~/do_heart/results/global_enrichment_allez/"
# set up directory (HPC)
enrich_dir  <- "/Users/c-alfond/do_heart/results/global_enrichment_allez/"
setwd(enrich_dir)

#########################################
#########################################

# Helper function to create the Jaccards distance df
makeJaccard <- function(data){
  
  print(data)
  # create list of gene lists
  parsed_genelist <- vector("list", length = nrow(data))
  for(n in 1:nrow(data)){
    parsed_genelist[[n]] <- strsplit(as.character(data$symbols[n]), ",")
  }
  
  # create empty df to fill
  Jaccard_df <- as.data.frame(matrix(data = NA, nrow = nrow(data), ncol = nrow(data)))
  rownames(Jaccard_df) <- colnames(Jaccard_df) <- data$X
  
  # All vs All comparison
  for(x in 1:nrow(data)){
    print(x)
    for(y in 1:nrow(data)){
      AnB <- length(intersect(parsed_genelist[[x]][[1]], parsed_genelist[[y]][[1]]))
      AuB <- length(union(parsed_genelist[[x]][[1]], parsed_genelist[[y]][[1]]))
      Jaccard_df[x,y] <- AnB / AuB
    }
  }
  return(Jaccard_df)
}


#########################################
#########################################

# get list of all files, calculate their Jaccard disance and make plots (in parallel)

file_list <- dir(pattern = "gsea_tbl_allez.csv")
types <- c("mrna_age", "mrna_sex", "protein_age", "protein_sex")

######

# Plot heatmaps (in parallel)
pdf("allez_categories_reduce_heatmaps.pdf")
# foreach(type = types) %dopar% {
for (type in types) {
  # Debugging
  file <- file_list[grep(type, file_list)]
  # output base file name
  basename <- gsub('.csv', '', file)
  
  print(basename)
  data <- read.csv(file)
  
  # Skip file if nothing is in it.
  if (nrow(data) == 0) {next}
  # Calculate Jaccard distance. 
  Jaccard_df <- makeJaccard(data)
  
  # plot heatmap
  heatmap.2(as.matrix(Jaccard_df), trace = "none", density.info = "none", 
            key=FALSE, cexRow = 1, cexCol = 1, mar=c(6,6), xlab = NULL, ylab = NULL)
  
  text(unit(0.60, "npc"), unit(0.77, "npc"), type)
  
  # # Plot dendogram
  # hc <- hclust(dist(as.matrix(Jaccard_df)))
  # plot(hc)
}
dev.off()


####################
####################
# Beyond this point, manually determine number of groups visually by plots for each file.
###########################
# 'my_cluster_func': Helper function to create clusters of categories for several files at a time.
# file_list is a list of files containing enriched category data that should be clustered.
# protein_ageN, protein_sexN, etc are numbers indicating how many clusters should be made out of each data file.
my_cluster_func <- function(file_list,  mrna_ageN, mrna_sexN, protein_ageN, protein_sexN){
  
  for(type in c("mrna_age", "mrna_sex", "protein_age", "protein_sex")){
    
    # get file and its data
    file <- file_list[grepl(type, file_list)]
    data <- read_csv(file)
    
    # Calculate Jaccard distance
    Jdist <- makeJaccard(data)
    hc <- hclust(dist(as.matrix(Jdist)))
    
    # Cluster according to file type. Each file type should can
    # get a different number of clusters according to the parameter 
    # values passed to this function
    
    if(type == "mrna_age"){
      if(is.na(mrna_ageN)){
        next
      }
      clust <- cutree(hc, k = mrna_ageN)
      
    } else if(type == "mrna_sex"){
      if(is.na(mrna_sexN)){
        next
      }
      clust <- cutree(hc, k = mrna_sexN)
      
    } else if (type == "protein_age"){
      if(is.na(protein_ageN)){
        next
      }
      clust <- cutree(hc, k = protein_ageN)
      
    } else if(type == "protein_sex"){
      if(is.na(protein_sexN)){
        next
      }
      clust <- cutree(hc, k = protein_sexN)
    }
    
    # Assign cluster group to each category.
    data$ClusterN <- clust
    # Order by group
    data <- arrange(data, ClusterN)
    
    # save data
    basename <- gsub(".csv", "", file)
    write.csv(data, file = paste0(basename, "_clustered.csv"), row.names = FALSE)
  }
}


file_list <- c("mrna_age_gsea_tbl_allez.csv",  "mrna_sex_gsea_tbl_allez.csv",
               "protein_age_gsea_tbl_allez.csv", "protein_sex_gsea_tbl_allez.csv")

# mrna_age has 28 groups
# mrna_sex has 18 groups
# protein_age has 4 groups
# protein_sex has 13 groups
file_list <- dir(pattern = "gsea_tbl_allez.csv")
my_cluster_func(file_list, mrna_ageN = 9, mrna_sexN = 4, protein_ageN = 4, protein_sexN = 4)

