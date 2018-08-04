#!/bin/bash
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=30
#PBS -N QTL_hotspots_job using data cmd
module load R/3.4.4
Rscript "/fastscratch/c-alfond/do_heart/scripts/pipeline/qtls/hotspots/QTL_hotspots.R"
