#!/bin/bash
#SBATCH --job-name=methylkit_diff
#SBATCH --account=project_2010912
#SBATCH --partition=small
#SBATCH --time=24:00:00          
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=64G                
#SBATCH --output=diff_meth_dss_ou%j.txt
#SBATCH --error=diff_meth_dss_er%j.txt

# Load r-env
module load r-env
echo "Job started on $(hostname)"
echo "Working directory: $(pwd)"
echo "CPUs: $SLURM_CPUS_PER_TASK"

cd /scratch/project_2010912/ezel/salmon
Rscript diff_meth_bf_correction.R