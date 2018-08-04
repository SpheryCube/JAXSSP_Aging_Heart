#!/bin/bash
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=38
#PBS -N qtl_allele_effects_stratified_job using data cmd
module load R/3.4.4
Rscript "/fastscratch/c-alfond/do_heart/scripts/pipeline/qtls/hotspots/qtl_allele_effects_stratified_v3.R"