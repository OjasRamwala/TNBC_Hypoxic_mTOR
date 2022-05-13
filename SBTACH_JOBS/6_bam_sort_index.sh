#!/bin/bash
#SBATCH --array=0-7
#SBATCH --time=10:00:00
#SBATCH --mem=64GB
#SBATCH --job-name=bam_sort_index
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END
#SBATCH --mail-user=smm924@nyu.edu
#SBATCH --output=bam_sort_index_%A_%a.out
#SBATCH --error=bam_sort_index_%A_%a.err


module purge
module load samtools/intel/1.14
SAM=($(ls *_aligned.sam))
SAM_FILES=${SAM[$SLURM_ARRAY_TASK_ID]}

samtools view --bam ${SAM_FILES} > ${SAM_FILES}.bam

BAM=($(ls *.bam))
BAM_FILES=${BAM[$SLURM_ARRAY_TASK_ID]}
samtools sort ${BAM_FILES} > ${BAM_FILES}.sorted.bam


SORTEDBAM=($(ls *.sorted.bam))
SORTEDBAM_FILES=${SORTEDBAM[$SLURM_ARRAY_TASK_ID]}

samtools index ${SORTEDBAM_FILES}
