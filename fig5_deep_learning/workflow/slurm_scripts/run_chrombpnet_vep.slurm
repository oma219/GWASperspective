#!/bin/bash

#SBATCH -J run_chrombpnet_vep
#SBATCH -o run_chrombpnet_vep.%A.stdout.txt
#SBATCH -e run_chrombpnet_vep.%A.stderr.txt
#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --time 00-16:00:00
#SBATCH --gres=gpu:1

SNP_GROUP_NAME="stingseq_nonhits"

REPO_DIR="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools"
OUTPUT="variant_scores/chrombpnet/${SNP_GROUP_NAME}_chrombpnet_vep_scores.csv"

printf "\n[log] running chrombpnet: generating ${OUTPUT}\n"
printf "[log] started running at $(date)\n\n"

cd $REPO_DIR
mamba run -n snakemake snakemake --use-conda \
                                 -c1 $OUTPUT

printf "[log] finished running chrombpnet at $(date)\n\n"