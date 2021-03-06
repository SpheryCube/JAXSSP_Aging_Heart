Jackson Labs Summer Student Program 2018
Daniel Alfonsetti, Massachusetts Institute of Technology
daniel.alfonsetti@gmail.com
alfonset@mit.edu
August 8, 2018

----------------
Project Overview
----------------

The goal of this project was to make a pipeline of statistical analyses that would identify how genetic variation affects the aging heart using transcriptome and proteome wide heart expression data. We had heart transcriptome (21101 measured transcripts) and proteome (4193 measured proteins) data for 189 DO mice, grouped by age of sacrifice (6, 12, and 18 months). Both sexes were included. QTL mapping, particularly interactive QTL mapping, was extensively used in this project. To get a better understanding of the project and its results, read 'final_report.pdf' or 'presentation.pptx' 

----------------
File structure
----------------
/data - Contains both original and cleaned data.

/scripts/pipeline - Contains scripts that perform a variety of analyses. Scripts that are part of a related analysis are grouped into subdirectories. Many R scripts have corresponding shell scripts so that they can also be run on the High Performance Computing Cluster (HPCC). Near the top of these scripts are areas where some code has to be commented or uncommented out depending on if you are using the HPC or not.

/scripts/exploratory - A couple scripts that dig deeper into the results outputted by the pipeline scripts.

/results - Contains all of the output of the scripts (besides the output of clean_data.R). Contains subdirectories that correspond to subdirectories in /scripts


----------------
Acknowledgements
----------------
Dan Gatti, Ph.D
Gary Churchill, Ph.D (Principal Investigator)
Susan McClatchy (Mentor)
Yuka Takemon
This Summer Student Fellowship was supported by Educational Gifts to The Jackson Laboratory


