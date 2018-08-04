#!/bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=25
#PBS -l mem=64GB
#PBS -N get_qtl_peaks_job using data cmd
module load R/3.4.4
Rscript "/fastscratch/c-alfond/do_heart/scripts/pipeline/qtls/scans/get_qtl_peaks.R"