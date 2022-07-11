#!/bin/bash
#SBATCH --job-name="rda rbt"
#SBATCH --cpus-per-task=40
#SBATCH --mem=256GB
#SBATCH --output=rda-model-job.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tseaborn@uidaho.edu
module load R/4.0.2
Rscript 05_prda.R