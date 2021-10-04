#!/bin/bash
#SBATCH --job-name="RBT RDA"
#SBATCH --cpus-per-task=47
#SBATCH --mem=1000GB
#SBATCH --output=long-job.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tseaborn@uidaho.edu
module load R/4.0.5
Rscript 04_env_prep.R
Rscript 05_rda.R