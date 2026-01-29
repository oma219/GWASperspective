#!/bin/bash

#SBATCH -J run_borzoi_prime
#SBATCH -o run_borzoi_prime.%A.stdout.txt
#SBATCH -e run_borzoi_prime.%A.stderr.txt
#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --time 02-00:00:00
#SBATCH --gres=gpu:1


REPO_DIR="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools"
OUTPUT="generate_all_possible_stingseq_hit_borzoi_prime_scores"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

printf "\n[log] started generating ${OUTPUT}\n"
printf "[log] started running at $(date)\n\n"

BORZOI_PRIME_MODEL="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools/data/borzoi_models/borzoi_prime/model0_best.h5"
BORZOI_PRIME_MODEL_PARAMS="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools/data/borzoi_models/borzoi_prime/params.json"

cd $REPO_DIR
snakemake --use-conda \
          --conda-frontend mamba \
            --config \
            BORZOI_MODEL=$BORZOI_PRIME_MODEL \
            BORZOI_MODEL_PARAMS=$BORZOI_PRIME_MODEL_PARAMS \
          -c1 $OUTPUT

printf "[log] finished running at $(date)\n\n"