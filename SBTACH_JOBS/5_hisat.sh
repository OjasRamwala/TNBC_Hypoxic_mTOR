#!/bin/bash
#SBATCH --array=0-7
#SBATCH --time=10:00:00
#SBATCH --mem=64GB
#SBATCH --job-name=HISAT
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END
#SBATCH --mail-user=or2072@nyu.edu
#SBATCH --output=HISAT_%A_%a.out
#SBATCH --error=HISAT_%A_%a.err

module purge
module load hisat2/2.2.1
FILES_1=($(ls *_1.fastq))
FILES_2=($(ls *_2.fastq))

INPUT_1=${FILES_1[$SLURM_ARRAY_TASK_ID]}
INPUT_2=${FILES_2[$SLURM_ARRAY_TASK_ID]}

echo $INPUT_1 $INPUT_2

hisat2 --threads 8 -x hg38 -1 ${INPUT_1} -2 ${INPUT_2} -S ${INPUT_1}_aligned.sam
