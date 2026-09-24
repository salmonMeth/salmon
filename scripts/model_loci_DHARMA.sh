#!/bin/bash
#SBATCH --job-name=glmm1000
#SBATCH --account=project_2010912
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=glmm_chr3_%j.out
#SBATCH --error=glmm_chr3_%j.err

module load r-env

Rscript /scratch/project_2010912/ezel/salmon/model_loci_overdisp_DHARMA.R