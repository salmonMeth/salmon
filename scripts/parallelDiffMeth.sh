#!/bin/bash
#SBATCH --job-name=methylkit_diff
#SBATCH --account=project_2010912
#SBATCH --partition=small
#SBATCH --time=24:00:00          
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=64G                
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt

# Load r-env
module load r-env

# Run the script
srun Rscript diffMeth.R

