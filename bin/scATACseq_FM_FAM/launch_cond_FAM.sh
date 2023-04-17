#!/bin/bash

#SBATCH --job-name=cond_markers
#SBATCH --partition long
#SBATCH --mem=40GB
#SBATCH --account=humancd34_diff_rna_atacseq
#SBATCH -c 6

module load r

Rscript /shared/projects/humancd34_diff_rna_atacseq/scATACseq/script/script_cond_FAM.R
