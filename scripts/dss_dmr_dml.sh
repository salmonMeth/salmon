#!/bin/bash
#SBATCH --job-name=DSS
#SBATCH --account=project_2010912
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=dss_%j.out
#SBATCH --error=dss_%j.err

module load r-env

Rscript /scratch/project_2010912/ezel/salmon/delete_later.R