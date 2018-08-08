/scripts/pipeline/



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

The order to run files in qtls/ is as follows:

1. scans/qtl.R
2. scans/get_qtl_peaks.R
3. hotspots/QTL_hotspots.R
4. hotspots/qtl_allele_effects_stratified.R
5. hotspots/mediation_analysis.R