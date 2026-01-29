#!/bin/bash

#SBATCH -J run_borzoi_prime_vep
#SBATCH -o run_borzoi_prime_vep.%A.stdout.txt
#SBATCH -e run_borzoi_prime_vep.%A.stderr.txt
#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --time 00-16:00:00
#SBATCH --gres=gpu:1

SNP_GROUP_NAME="stingseq_nonhits"

REPO_DIR="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools"
OUTPUT="variant_scores/borzoi_prime/${SNP_GROUP_NAME}_borzoi_prime_vep_scores.csv"

BORZOI_PRIME_MODEL="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools/data/borzoi_models/borzoi_prime/model0_best.h5"
BORZOI_PRIME_MODEL_PARAMS="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools/data/borzoi_models/borzoi_prime/params.json"

printf "\n[log] running borzoi prime: generating ${OUTPUT}\n"
printf "[log] started running at $(date)\n\n"

cd $REPO_DIR
mamba run -n snakemake snakemake --use-conda \
                                 --conda-frontend mamba \
                                 --config \
                                 BORZOI_MODEL=$BORZOI_PRIME_MODEL \
                                 BORZOI_MODEL_PARAMS=$BORZOI_PRIME_MODEL_PARAMS \
                                 -c1 $OUTPUT

printf "[log] finished running borzoi at $(date)\n\n"