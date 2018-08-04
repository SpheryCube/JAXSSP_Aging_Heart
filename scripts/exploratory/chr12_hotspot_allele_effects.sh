#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=8
#PBS -N qtl_chr12_interaction_genes_explore_job using data cmd
module load R/3.4.4
Rscript "/fastscratch/c-alfond/do_heart/scripts/exploratory/chr12_hotspot_allele_effects.R"


