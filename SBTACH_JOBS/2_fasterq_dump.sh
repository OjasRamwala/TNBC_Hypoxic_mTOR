#!/bin/bash
#SBATCH --cpus-per-task=15
#SBATCH --time=5:00:00
#SBATCH --mem=64GB
#SBATCH --mail-user=or2072@nyu.edu
#SBATCH --output=prefetch

module purge
module load sra-tools/2.10.9
fasterq-dump all_sra_files
