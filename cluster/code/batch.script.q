#!/bin/bash
#
#SBATCH -o Rtest.out  # will be empty if there's no output
#SBATCH -e Rtest.err        # useful for troubleshooting, will be empty if there are no errors
#SBATCH -N 1
#SBATCH -J "p0W0X0Y0Z0"
#SBATCH --array=1-19
#SBATCH -t 0-12:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ifmiller@princeton.edu


INDEX=$(( $SLURM_ARRAY_TASK_ID + 0 ))
cd ~/immuno.var.vir.evo/p0W0X0Y0Z0/dir.$INDEX
srun R CMD BATCH analysis.R
