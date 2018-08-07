/results/hotspots/

-----------------
File Structure
-----------------
This folder contains the output from the scripts in /scripts/pipeline/qtls/hotspots/. There is a subfolder for 6 of the 10 peaks tables in results/ (all the tables except the tables for the full models). These subfolders are...

mrna_age_hotspots_output/ contains output for the pipeline analysis for age-interactive eQTLs 
mrna_none_hotspots_output/ contains output for the pipeline analysis for purely additive eQTLs
mrna_sex_hotspots_output/ contains output for the pipeline analysis for sex-interactive eQTLs

protein_*/ directories have the same meaning, but for pQTLs.


Additionally, there is the...

preliminary_hotspot_results/ subfolder which contains QTL maps and QTL density histograms for each of the 6 scan types being analyzed.

--------------------------
Inside a hotspot folder...
--------------------------

Inside each of these folders are files with outputs that analyze QTL hotspots for that folder's scan type. Some folders/scan types have multiple hotspots. I will use the files in the "protein_age_hotspots_output/" folder as an example to explain here. Each folder has similiar files.

Note: Only trans-hotspots were analyzed (It's not really possible to have a cis hotspot anyways (?))

protein_age_trans_chr3_hotspot_genes.csv - A file containing gene information and QTL information for all the genes that have QTL in the hotspot.

protein_age_trans_chr3_hotspot_go_gsea_cluster.csv / protein_age_trans_chr3_hotspot_kegg_gsea_cluster.csv - Files containing GO and KEGG categories that are enriched with genes that have QTL in the hotspot.

protein_age_trans_hotspot_pcs.pdf - A file displaying formation about the first 3 principal components generated from the expression data of genes that have QTL in the hotspot (a seperate page is made for each hotspot if there are multiple hotspots to be analyzed)

protein_age_trans_qtl_gene_hotspot_mrna_correlations.pdf - A file displaying correlation heatmaps for the mRNA expression of each outcome gene in the QTL.

protein_age_trans_qtl_gene_hotspot_protein_correlations.pdf - The same as above, but now using protein expression correlations.

protein_age_trans_pqtl_denisty_hist.pdf - An image showing where QTLs stack up along genome.

protein_age_QTLmap.pdf - An image associating a gene location's to that genes pQTL.

allele_effect_outputs/ - A folder containing an in-depth analysis of each gene that has a QTL in a hotspot.


