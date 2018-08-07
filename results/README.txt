/results

This folder contains outputs from the scripts in /scripts. In particular, it contains several subdirectories corresponding to subdirectories in in /scripts. QTL peak table files are in this current directory.

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