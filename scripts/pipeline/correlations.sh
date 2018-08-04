#!/bin/bash
#
#PBS -N correlations_job using data cmd
#PBS -q batch
#PBS -m bae
#PBS -l nodes=1:ppn=4

module load R/3.4.4
Rscript "/fastscratch/c-alfond/do_heart/scripts/pipeline/correlations.R"