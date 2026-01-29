#!/bin/bash

#SBATCH -J run_ag_on_stingseq_loci
#SBATCH -o run_ag_on_stingseq_loci.%A.stdout.txt
#SBATCH -e run_ag_on_stingseq_loci.%A.stderr.txt
#SBATCH -p bigmem
#SBATCH -c 4
#SBATCH --mem=100G
#SBATCH --time 00-16:00:00
#SBATCH --array=1-139%5

##############################################
# step 1: get input sumstats file to process
##############################################

WORK_DIR="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools/results"

INPUT_LIST_DIR="analyses/stingseq_datasets/all_hit_loci/metadata"
OUTPUT_DIR="analyses/stingseq_datasets/all_hit_loci/variant_scores"

FILE_LIST="${WORK_DIR}/${INPUT_LIST_DIR}/validated_sumstats_list.txt"

# important: I filter the filelist to relevant sumstats
input_sumstats_file=$(grep "\-ukbb-" "$FILE_LIST" | awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n {print; exit}')

if [[ -z "$input_sumstats_file" ]]; then
    echo "No line for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

filename=$(basename "$input_sumstats_file")
sumstats_id="${filename%.hg38.csv}"

OUTPUT="${OUTPUT_DIR}/${sumstats_id}.alphagenome_vep_scores.csv"

##############################################
# step 2: run alphagenome scoring
##############################################
REPO_DIR="/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/run_deep_learning_tools"
cd $REPO_DIR

source ~/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --use-conda \
          --conda-frontend mamba \
          -c1 $OUTPUT
