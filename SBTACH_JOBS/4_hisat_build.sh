#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=128GB
#SBATCH --job-name=HISAT
#SBATCH --mail-type=END
#SBATCH --mail-user=or2072@nyu.edu
#SBATCH --output=HISAT_BUILD.out

module purge
module load hisat2/2.2.1
hisat2-build hg38.fa hg38
