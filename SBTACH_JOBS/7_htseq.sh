#!/bin/bash
#SBATCH --array=0-7
#SBATCH --time=10:00:00
#SBATCH --mem=64GB
#SBATCH --job-name=HTSeq
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END
#SBATCH --mail-user=or2072@nyu.edu
#SBATCH --output=HTSeq_%A_%a.out
#SBATCH --error=HTSeq_%A_%a.err

module purge
module load htseq/0.13.5

BAM=($(ls *.sorted.bam))
BAM_FILES=${BAM[$SLURM_ARRAY_TASK_ID]}

htseq-count ${BAM_FILES} hg38.ensGene.gtf > ${BAM_FILES}.count
