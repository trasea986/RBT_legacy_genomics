#!/bin/bash
#SBATCH --job-name="RBT LFMM"
#SBATCH --cpus-per-task=24
#SBATCH --mem=192GB
#SBATCH --output=lfmm-model-job.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tseaborn@uidaho.edu
module load R/4.0.2
Rscript 08_lfmm.R