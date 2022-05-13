#!/bin/bash
#SBATCH --cpus-per-task=15
#SBATCH --time=5:00:00
#SBATCH --mem=64GB
#SBATCH --mail-user=or2072@nyu.edu
#SBATCH --output=fasterq_dump

module purge
module load sra-tools/2.10.9
fasterq-dump *.sra

# fasterq-dump SRR6068874
# fasterq-dump SRR6068875
# fasterq-dump SRR6068876
# fasterq-dump SRR6068877
# fasterq-dump SRR6068882
# fasterq-dump SRR6068883
# fasterq-dump SRR6068884
# fasterq-dump SRR6068885
# fasterq-dump SRR6068890
# fasterq-dump SRR6068891
# fasterq-dump SRR6068892
# fasterq-dump SRR6068893
# fasterq-dump SRR6068898
# fasterq-dump SRR6068899
# fasterq-dump SRR6068900
# fasterq-dump SRR6068901
