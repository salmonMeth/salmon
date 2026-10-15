#!/bin/bash
#SBATCH --account=project_2010912
#SBATCH --job-name=simulate_loci_model
#SBATCH --output=logs/sim_model_loci_%j.out
#SBATCH --error=logs/sim_model_loci_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

module load r-env

# Run the R script
Rscript /scratch/project_2010912/ezel/salmon/sim_loci_model.R
