#!/bin/bash
#SBATCH --job-name="RBT LFMM"
#SBATCH --cpus-per-task=12
#SBATCH --mem=192GB
#SBATCH --output=lfmm-sim-job.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tseaborn@uidaho.edu
module load R/4.0.2
Rscript 07_lfmm_sim_ind.R