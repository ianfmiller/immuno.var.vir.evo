#!/bin/bash
#
#SBATCH -o data.out  # will be empty if there's no output
#SBATCH -e data.err        # useful for troubleshooting, will be empty if there are no errors
#SBATCH -N 1
#SBATCH -J "p0W0X0Y0Z0"
#SBATCH --array=1-1
#SBATCH -t 0-00:05:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ifmiller@princeton.edu

cd ~/immuno.var.vir.evo/p0W0X0Y0Z0
srun R CMD BATCH data.compilation.R
