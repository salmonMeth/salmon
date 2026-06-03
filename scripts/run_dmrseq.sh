#!/bin/bash
#SBATCH --job-name=dmrseq
#SBATCH --account=project_2010912
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=dmrseq.out
#SBATCH --error=dmrseq.err

module load r-env

Rscript run_dmrseq.R