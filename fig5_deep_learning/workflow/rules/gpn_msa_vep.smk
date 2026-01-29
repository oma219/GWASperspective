######################################################
# Name: gpn_msa_vep.smk
# Description: score variants using gpn msa
# Date: January 3rd, 2026
######################################################

rule run_gpn_msa_vep:
    input:
        variant_list="preprocessed_snp_lists/{group_name}_snp_list.csv"
    output:
        scores="variant_scores/gpn_msa/{group_name}_gpn_msa_vep_scores.csv"
    shell:
        """
        ########################################################
        # section 0: load tabix to use for parsing scores file
        ########################################################
        module load all/tabixpp/1.1.2-GCC-12.3.0

        ###############################################
        # section 1: iterate through variants and pull 
        # out score from the precomputed scores file
        ###############################################
        printf "chrom,pos,ref,alt,tool,score\n" > {output.scores}

        num_variants_processed=0
        while IFS=, read -r chrom pos ref alt; do
            # skip the header line
            if [[ "$chrom" == "chrom" ]]; then
                continue
            fi

            # check if all 4 fields are non-empty
            if [[ -z "$chrom" || -z "$pos" || -z "$ref" || -z "$alt" ]]; then
                echo "[error] line does not have 4 fields: $chrom,$pos,$ref,$alt" 1>&2
                exit 1
            fi

            # pull out the score
            score=$(tabix {gpn_msa_hg38_scores} "${{chrom}}:${{pos}}-${{pos}}" | awk -F"\t" -v a="$alt" '$4 == a {{print ($5 < 0 ? -$5 : $5)}}')
            printf "$chrom,$pos,$ref,$alt,gpn_msa,$score\n" >> {output.scores}
            num_variants_processed=$((num_variants_processed + 1))

        done < {input.variant_list}
        
        printf "[log] finished pulling out the scores for each variant (n = $num_variants_processed variants)\n\n"
        """
