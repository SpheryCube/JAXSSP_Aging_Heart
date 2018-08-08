/scripts/qtls/scans/

qtl.R - Runs QTLs using additive and full models (additive covs + int cover), saves the results, and then subtracts them in order to get interaction effect. Permutation testing can be run in this file too, but currently we are just using a flat threshold of 6 since we are going to be looking for hotspots.

get_qtl_peaks.R - Pulls out peaks for each type of scan. Currently using a LOD threshold of 6.



Each of these files have shell scripts so that they can be run on the HPC. Some changes need to be made at the top of each R file (by comment or uncommenting) in order to be usable on the HPC.