/scripts/pipeline/

Scripts in this folder 

1. clean_data.R

There are two seperate pipelines. The first one does correlations, ANOVAs, and gene set enrichment:

2. anova/anova.R
3. anova/mrna_corrs_v_protein_corrs.R
5. mRNA_protein_cors.R
6. mRNA_v_protein_by_age.R
7. global_enrichment_allez.R
8. allez_category_reduce.R

The second one analyses additive and interactive QTL mapping, QTL hotspots, and allele/genotype effects.
It is broken down into two folders: scans/ and hotspots/

scans/ -  A folder containing files to run additivie and interactive QTL scans and then extract peak tables.
hotspots/ - A folder containing scripts to find hotspots regions in each type of QTL scan and then analyze the genes that have QTLs in those hotspots.
