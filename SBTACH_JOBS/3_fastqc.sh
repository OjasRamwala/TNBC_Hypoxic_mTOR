#!/bin/bash
#SBATCH --array=0-3
#SBATCH --time=10:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=FASTQC
#SBATCH -n 4
#SBATCH --mail-type=END
#SBATCH --mail-user=or2072@nyu.edu
#SBATCH --output=FASTQC_%A_%a.out
#SBATCH --error=FASTQC_%A_%a.err

module purge #removing all the previously loaded modules
module load fastqc/0.11.9 #loading the corresponding fastqc module

# capture the output of a command line and store it in a variable
FILES=($(ls *.fastq))

# this is going to assign the variables to file names
INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}

fastqc $INPUT
