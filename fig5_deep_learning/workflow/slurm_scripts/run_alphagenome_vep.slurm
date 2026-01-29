#!/bin/bash

#SBATCH -J run_ag_vep
#SBATCH -o run_ag_vep.%A.stdout.txt
#SBATCH -e run_ag_vep.%A.stderr.txt
#SBATCH -p bigmem
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --time 00-08:00:00

SNP_GROUP_NAME="stingseq_nonhits"

REPO_DIR="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools"
OUTPUT="variant_scores/alphagenome/${SNP_GROUP_NAME}_alphagenome_vep_scores.csv"

printf "\n[log] running alphagenome: generating ${OUTPUT}\n"
printf "[log] started running at $(date)\n\n"

cd $REPO_DIR
mamba run -n snakemake snakemake --use-conda \
                                 --conda-frontend mamba \
                                 -c1 $OUTPUT

printf "[log] finished running alphagenome at $(date)\n\n"