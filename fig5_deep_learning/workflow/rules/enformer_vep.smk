######################################################
# Name: enformer_vep.smk
# Description: score variants using enformer
# Date: December 26th, 2025
######################################################

"""
example of running this rule:

    snakemake --use-conda \
              -c1 variant_scores/enformer/{group_name}_enformer_vep_scores.csv

notes:
    - notice the modules that are loaded prior to running the script,
      this code is specific to the NYGC cluster. 
"""

######################################################
# section 0: basic version used for simple cases
######################################################

rule run_enformer_vep:
    input:
        variant_list="preprocessed_snp_lists/{group_name}_snp_list.csv"
    output:
        scores="variant_scores/enformer/{group_name}_enformer_vep_scores.csv"
    conda:
        "../../envs/enformer_py38.yaml"
    shell:
        """
        module load cuda/11.7.0 
        module load cudnn/8.5.0.96-CUDA-11.7.0

        python3 {repo_dir}/src/enformer_vep.py --variant-list {input.variant_list} \
                                               --genome {hg38_genome} \
                                               --output {output.scores}
        """

######################################################
# section 1: specific rule for running scoring on
#            sting-seq hit loci sumstats
###################################################### 

rule run_enformer_vep_stingseq_hit_loci:
    input:
        variant_list="analyses/stingseq_datasets/all_gwas_loci/sumstats/{loci_type}/variant{id}-{trait}-{consort}-variants.hg38.csv"
    output:
        scores="analyses/stingseq_datasets/all_gwas_loci/variant_scores/enformer/{loci_type}/variant{id}-{trait}-{consort}-variants.enformer_vep_scores.csv"
    conda:
        "../../envs/enformer_py38.yaml"
    shell:
        """
        module load cuda/11.7.0 
        module load cudnn/8.5.0.96-CUDA-11.7.0

        python3 {repo_dir}/src/enformer_vep.py --variant-list {input.variant_list} \
                                               --genome {hg38_genome} \
                                               --output {output.scores}
        """