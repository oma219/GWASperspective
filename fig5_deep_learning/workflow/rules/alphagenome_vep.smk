######################################################
# Name: alphagenome_vep.smk
# Description: score variants using alphagenome
# Date: December 17th, 2025
######################################################

"""
example of running this rule:

    snakemake --use-conda \
              -c1 variant_scores/alphagenome/{group_name}_alphagenome_vep_scores.csv

notes:
    - this command will use the alphagenome api key
      defined in the config file to change it, you can specify
      it on the command line.
"""

######################################################
# section 0: basic version used for simple cases
######################################################

rule run_alphagenome_vep:
    input:
        variant_list="preprocessed_snp_lists/{group_name}_snp_list.csv"
    output:
        scores="variant_scores/alphagenome/{group_name}_alphagenome_vep_scores.csv"
    conda:
        "../../envs/alphagenome_py311.yaml"
    shell:
        """
        python3 {repo_dir}/src/alphagenome_vep.py --variant-list {input.variant_list} \
                                                  --api-key {alphagenome_api_key} \
                                                  --seq-length 1MB \
                                                  --output {output.scores}
        """

######################################################
# section 1: specific rule for running scoring on
#            sting-seq hit loci sumstats
###################################################### 

rule run_alphagenome_vep_stingseq_hit_loci:
    input:
        variant_list="analyses/stingseq_datasets/all_gwas_loci/sumstats/{loci_type}/variant{id}-{trait}-{consort}-variants.hg38.csv"
    output:
        scores="analyses/stingseq_datasets/all_gwas_loci/variant_scores/alphagenome/{loci_type}/variant{id}-{trait}-{consort}-variants.alphagenome_vep_scores.csv"
    conda:
        "../../envs/alphagenome_py311.yaml"
    shell:
        """
        python3 {repo_dir}/src/alphagenome_vep.py --variant-list {input.variant_list} \
                                                  --api-key {alphagenome_api_key} \
                                                  --seq-length 1MB \
                                                  --output {output.scores}
        """
