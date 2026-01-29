######################################################
# Name: chrombpnet_vep.smk
# Description: score variants using chrombpnet
# Date: January 3rd, 2026
######################################################

rule run_chrombpnet_vep:
    input:
        variant_list="preprocessed_snp_lists/{group_name}_snp_list.csv"
    output:
        fold_scores="variant_scores/chrombpnet/{group_name}_chrombpnet_vep_scores.fold_{fold_i}.variant_scores.tsv"
    conda:
        "../../envs/chrombpnet_py38.yaml"
    shell:
        """
        #######################################################
        # section 0: create a temporary file for reorganizing 
        #            the variant list to match scheme expected 
        #            by chrombpnet
        #######################################################
        tmpfile="$(mktemp)"
        trap 'rm -f "$tmpfile"' EXIT

        awk -F, 'NR>1 {{print $1"\t"$2"\t"$3"\t"$4"\t"$1"_"$2"_"$3"_"$4}}' {input.variant_list} > $tmpfile
        cat $tmpfile

        #######################################################
        # section 1: run chrombpnet scoring for this fold
        #######################################################
        module load cuda/11.7.0
        module load cudnn/8.5.0.96-CUDA-11.7.0

        bias_model="{data_dir}/chrombpnet_models/fold_{wildcards.fold_i}/model.bias_scaled.fold_{wildcards.fold_i}.ENCSR868FGK.h5"
        nobias_model="{data_dir}/chrombpnet_models/fold_{wildcards.fold_i}/model.chrombpnet_nobias.fold_{wildcards.fold_i}.ENCSR868FGK.h5"
        output_prefix="variant_scores/chrombpnet/{wildcards.group_name}_chrombpnet_vep_scores.fold_{wildcards.fold_i}"

        python3 ~/variant-scorer/src/variant_scoring.py -l $tmpfile\
                                                        -g {hg38_genome} \
                                                        -s {hg38_chrom_sizes} \
                                                        -m $nobias_model \
                                                        -b $bias_model \
                                                        -sc chrombpnet \
                                                        -o $output_prefix
        """

rule run_chrombpnet_ensemble_vep_averaging:
    input:
        fold_scores=lambda wc: expand(
            "variant_scores/chrombpnet/{group_name}_chrombpnet_vep_scores.fold_{fold_i}.variant_scores.tsv",
            group_name=wc.group_name,
            fold_i=range(chrombpnet_num_folds),
        )
    output:
        avg_scores="variant_scores/chrombpnet/{group_name}_chrombpnet_vep_scores.all_folds.mean.variant_scores.tsv"
    conda:
        "../../envs/chrombpnet_py38.yaml"
    shell:
        """
        output_prefix="variant_scores/chrombpnet/{wildcards.group_name}_chrombpnet_vep_scores.all_folds"

        score_list=""
        score_dir="variant_scores/chrombpnet/"
        for fold_i in $(seq 0 $(({chrombpnet_num_folds} - 1))); do
            curr_score_file="{wildcards.group_name}_chrombpnet_vep_scores.fold_${{fold_i}}.variant_scores.tsv"
            score_list="$score_list $curr_score_file"
        done

        python3 ~/variant-scorer/src/variant_summary_across_folds.py \
                                        --score_dir $score_dir \
                                        --score_list $score_list \
                                        --out_prefix $output_prefix \
                                        --schema chrombpnet
        """

rule parse_chrombpnet_all_folds_output:
    input:
        avg_scores="variant_scores/chrombpnet/{group_name}_chrombpnet_vep_scores.all_folds.mean.variant_scores.tsv"
    output:
        scores="variant_scores/chrombpnet/{group_name}_chrombpnet_vep_scores.csv"
    run:
        import pandas as pd 

        scores_df = pd.read_csv(input.avg_scores, sep='\t', dtype=str)

        with open(output.scores, "w") as out_fd:
            out_fd.write("chrom,pos,ref,alt,tool,score\n")

            for index, row in scores_df.iterrows():
                chrom = row['chr'].replace('chr', '')
                pos = row['pos']
                ref = row['allele1']
                alt = row['allele2']
                tool = 'chrombpnet'
                score = row['abs_logfc.mean']

                out_fd.write(f"{chrom},{pos},{ref},{alt},{tool},{score}\n")