#!/bin/bash
#SBATCH --account=project_2010912
#SBATCH --job-name=simulate_methylation
#SBATCH --output=logs/simulate_%j.out
#SBATCH --error=logs/simulate_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

module load r-env

# Run the R script
Rscript /scratch/project_2010912/ezel/salmon/sim_loci_meth.R

