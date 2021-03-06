/results

This folder contains outputs from the scripts in /scripts. In particular, it contains several subdirectories corresponding to subdirectories in in /scripts. QTL peak table files are in this current directory.
------------------------------
File structure
------------------------------
anova/ - results from ANOVA scripts and the gene set enrichment analysis done on each quadrant of the mRNA-Age vs protein-Age correlations plot

global_enrichment_allez/ - Enrichment analysis results after performing differential expression on age and differential expression on sex

hotspots/ - Results from the hotspots output. Has many subdirectories

mrna_protein_cors/ - Correlation outputs between mrna and their respective protein. Gives raw, age-adjusted, and sex-adjusted results.

mrna_protein-cors_by_age_sex/ - We looked to see if the correlations between mrna and their respective proteins changed with age. We found they did (they decreased on average), but the way in which they decreased was depenendent on sex. 


------------------------------
Peak file title interpretation:
------------------------------

peaks_mrna_none_thresh_6.csv - eQTL scans were only done with additive covariates. Only eQTL peaks with a LOD of 6 or above were recorded.

peaks_mrna_age_full_thresh_6.csv eQTL scans were done with additive covariates plus age as the interactive covariate

peaks_mrna_age_thresh_6.csv - QTL-age interaction effect peaks. Gotten by subtracting the purely additive QTL model output from the full QTL age interactive model output (before taking the threshold). Only interaction effect peaks with a LOD of 6 or above were taken.

peaks_mrna_sex_full_thresh_6.csv - eQTL scans done with additive covariates plus sex as the interactive covariate.

peaks_mrna_sex_thresh_6.csv - QTL-sex interaction effect peaks. Gotten by subtracting the purely additive eQTL model output from the full QTL sex interactive model output (before taking the threshold). Only interaction effect peaks with a LOD of 6 or above were taken.


peaks_protein_... files have the same meaning except they are pQTLs instead of eQTLs.

Note: Additive covariates for each scan were sex and age.
------------------------------