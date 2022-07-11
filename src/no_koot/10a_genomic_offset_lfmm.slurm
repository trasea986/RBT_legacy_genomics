#!/bin/bash
#SBATCH --job-name="RBT LFMM Genomic Offset"
#SBATCH --cpus-per-task=24
#SBATCH --mem=192GB
#SBATCH --output=lfmm-genomic-offset.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tseaborn@uidaho.edu
module load R/4.0.2
Rscript 10a_genomic_offset_cluster_model.R