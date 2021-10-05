# Host heterogeneity and the evolution of quantitative pathogen virulence
## To run all analyses on a cluster with SLURM based job scheduling:
### 1) Copy the contents of the 'cluster' file to a directory titled 'immuno.var.vir.evo' in your home directory.
### 2) Enter the command 'Rscript run_analysis.R'
#### For each combination of w,x,y, and z values, this script
##### 1) Creates a new subdirectory
##### 2) Creates a copy of the analysis script for each immunity distribution shape parameter (theta), set to run for the w,x,y,z, and theta parameters of interest
##### 3) Creates a copy of the data compilation script for the w,x,y, and z parameters of interest
##### 4) Writes batch scripts for running the analyses (using "analysis.R") and compiling output (using "data.compilation.R")
##### 5) Executes these batch scripts.
### For a combination of w,x,y,z and theta parameters, the analysis script creates the epidemiological model using the "writeSIResscDD.R" script, conducts the adaptive dynamics analysis, and writes the result to a .RDS file.
### For a combination of w,x,y,z parameters, the data compliation script compiles all of the data (across the various values of theta) and writes the output as a .csv file
## To run an analysis on a personal machine, use the "analysis.R" script.
