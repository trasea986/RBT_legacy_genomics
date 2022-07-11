#!/bin/bash
#SBATCH --job-name="prda rbt"
#SBATCH --cpus-per-task=24
#SBATCH --mem=192GB
#SBATCH --output=prda-model-job.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tseaborn@uidaho.edu
module load R/4.0.2
Rscript 06_rda_viz.R