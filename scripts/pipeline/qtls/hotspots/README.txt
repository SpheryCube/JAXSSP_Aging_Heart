/do_heart/scripts/pipeline/qtls/hotspots/


QTL_hotspots.R - A two part script for finding QTL hotspots and then analyzing each one. Halfway through the script requires the user to look at the outputs from the first half and then decide what stacks of QTLs should actually be considered hotspots and further analyzed by specificy a QTL count threshold.

qtl_allele_effects_stratified.R - A script for analyzing the genes in each hotspot. Each gene will get its own output "story" file with a bunch of graphs. If the hotspot is for interactive QTL, the script will look at allele and genotype effects stratified by the levels of the interactive covariate.

mediation_analysis.R - Performs mediation analysis for interactive QTL scans using genes with QTL in the hotspots as targets and genes underneath/around the QTL location as candidate mediators. Both mrna and protein candidate mediators are used. Output files with _r10 at the end means candidate mediators were taken from a range of 10 Mbp around the target gene, while _r5 means candidate mediators were taken from a range of 5 Mbp around the target gene.

