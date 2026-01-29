######################################################
# Name: borzoi_prime_vep.smk
# Description: score variants using borzoi
# Date: January 2nd, 2025
######################################################

"""
example of running this rule:

    snakemake --use-conda \
              --config \
              BORZOI_MODEL=$BORZOI_PRIME_MODEL \
              BORZOI_MODEL_PARAMS=$BORZOI_PRIME_MODEL_PARAMS \
              -c1 variant_scores/borzoi_prime/{group_name}_borzoi_prime_vep_scores.csv

notes:
    - this command will use the borzoi prime mdoel and params
      but it is necessary to set them on the command line using
      the config parameter so it gets the prime model and not 
      the original model.
    - notice the environment variables that are modified prior to
      running the script, this was necessary to get approriate visibility
      for the cuda libraries to run the model.
"""

######################################################
# section 0: basic version used for simple cases
######################################################

rule run_borzoi_prime_vep:
    input:
        variant_list="preprocessed_snp_lists/{group_name}_snp_list.csv"
    output:
        scores="variant_scores/borzoi_prime/{group_name}_borzoi_prime_vep_scores.csv"
    conda:
        "../../envs/borzoi_py310.yaml"
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$CONDA_PREFIX/lib64:$LD_LIBRARY_PATH
        export CUDA_VISIBLE_DEVICES=0
        export XLA_FLAGS="--xla_gpu_cuda_data_dir=$CONDA_PREFIX"
        
        python3 {repo_dir}/src/borzoi_vep.py --variant-list {input.variant_list} \
                                             --model-params {borzoi_model_params} \
                                             --model {borzoi_model} \
                                             --genome {hg38_genome} \
                                             --tool-name borzoi_prime \
                                             --output {output.scores}
        """

######################################################
# section 1: specific rule for running scoring on
#            sting-seq hit loci sumstats
###################################################### 

rule run_borzoi_prime_vep_stingseq_hit_loci:
    input:
        variant_list="analyses/stingseq_datasets/all_gwas_loci/sumstats/{loci_type}/variant{id}-{trait}-{consort}-variants.hg38.csv"
    output:
        scores="analyses/stingseq_datasets/all_gwas_loci/variant_scores/borzoi_prime/{loci_type}/variant{id}-{trait}-{consort}-variants.borzoi_prime_vep_scores.csv"
    conda:
        "../../envs/borzoi_py310.yaml"
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$CONDA_PREFIX/lib64:$LD_LIBRARY_PATH
        export CUDA_VISIBLE_DEVICES=0
        export XLA_FLAGS="--xla_gpu_cuda_data_dir=$CONDA_PREFIX"
        
        python3 {repo_dir}/src/borzoi_vep.py --variant-list {input.variant_list} \
                                             --model-params {borzoi_model_params} \
                                             --model {borzoi_model} \
                                             --genome {hg38_genome} \
                                             --tool-name borzoi_prime \
                                             --output {output.scores}
        """