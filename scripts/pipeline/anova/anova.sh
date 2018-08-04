#!/bin/bash
#
#PBS -N ANOVAs_job using data cmd
#PBS -q batch
#PBS -m bae
#PBS -l nodes=1:ppn=4

module load R/3.4.4
Rscript "/fastscratch/c-alfond/do_heart/scripts/pipeline/anova/anova.R"