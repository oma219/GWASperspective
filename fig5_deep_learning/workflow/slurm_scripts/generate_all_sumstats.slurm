#!/bin/bash

#SBATCH -J generate_all_nonhit_sumstats
#SBATCH -o generate_all_nonhit_sumstats.%A.stdout.txt
#SBATCH -e generate_all_nonhit_sumstats.%A.stderr.txt
#SBATCH -p bigmem
#SBATCH -c 48
#SBATCH --mem=500G
#SBATCH --time 00-12:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=oahmed@nygenome.org

#################################################
# Note to self: this step took about ~50 minutes
#               to finish on bigmem node using
#               24 threads
#################################################
REPO_DIR="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools"
OUTPUT="generate_all_possible_stingseq_nonhit_hg19_sumstats"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

printf "\n[log] started generating ${OUTPUT}\n"
printf "[log] started running at $(date)\n\n"

cd $REPO_DIR
snakemake -c48 $OUTPUT

printf "[log] finished running at $(date)\n\n"