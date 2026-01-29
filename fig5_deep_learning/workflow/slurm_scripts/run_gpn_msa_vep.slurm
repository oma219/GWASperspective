#!/bin/bash

#SBATCH -J run_gpn_msa_vep
#SBATCH -o run_gpn_msa_vep.%A.stdout.txt
#SBATCH -e run_gpn_msa_vep.%A.stderr.txt
#SBATCH -p bigmem
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --time 00-16:00:00

SNP_GROUP_NAME="stingseq_nonhits"

REPO_DIR="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools"
OUTPUT="variant_scores/gpn_msa/${SNP_GROUP_NAME}_gpn_msa_vep_scores.csv"

printf "\n[log] running gpn msa: generating ${OUTPUT}\n"
printf "[log] started running at $(date)\n\n"

cd $REPO_DIR
mamba run -n snakemake snakemake -c1 $OUTPUT

printf "[log] finished running gpn msa at $(date)\n\n"