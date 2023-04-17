#!/bin/bash

#SBATCH --job-name=clust_markers
#SBATCH --partition long
#SBATCH --mem-per-cpu=40GB
#SBATCH --account=humancd34_diff_rna_atacseq
#SBATCH --ntasks=1

#SBATCH --array=0-15

module load r

srun Rscript /shared/projects/humancd34_diff_rna_atacseq/scATACseq/script/script_clust_FM.R $SLURM_ARRAY_TASK_ID
