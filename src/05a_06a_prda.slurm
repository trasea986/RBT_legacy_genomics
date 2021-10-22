#!/bin/bash
#SBATCH --job-name="prda rbt"
#SBATCH --cpus-per-task=32
#SBATCH --mem=192GB
#SBATCH --output=prda-model-job.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tseaborn@uidaho.edu
module load R/4.0.2
Rscript 05a_prda.R
Rscript 06a_prda_viz.R