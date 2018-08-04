#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=35
#PBS -l mem=64GB
#PBS -N mediation_analysis_job using data cmd
module load R/3.4.4
Rscript "/fastscratch/c-alfond/do_heart/scripts/pipeline/qtls/hotspots/mediation_analysis.R"
