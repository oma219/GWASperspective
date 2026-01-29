######################################################
# Name: extract_specific_stingseq_gwas_loci.smk
# Description: this file will extract the p-values,
#              pips, deep-learning percentiles for
#              a specific locus. 
#
#              important: this rule is not generalized 
#                         to work for any locus, just
#                         specific ones where the data
#                         has been gathered
#
# Date: january 16th, 2025
######################################################

rule build_masterfile_for_example_locus_1:
    output:
        masterfile="analyses/stingseq_datasets/example_loci/locus_1/example_locus_1_data.csv"
    run:
        import os
        import sys
        import numpy as np
        import pandas as pd
        from liftover import get_lifter
        
        converter = get_lifter('hg19', 'hg38', one_based=True)

        ##################################################
        # step 0: define helper methods
        ##################################################
        def lift_from_hg19_to_hg38(hg19_chrom, hg19_pos):
            pos_lifted = converter[f"chr{hg19_chrom}"][int(hg19_pos)]
            if len(pos_lifted) == 0:
                return False, "", ""
            else:
                hg38_chrom = pos_lifted[0][0].replace("chr", "")
                hg38_pos = int(pos_lifted[0][1])
                if hg38_chrom != hg19_chrom:
                    return False, "", ""

                return True, hg38_chrom, hg38_pos

        ############################################################
        # step 1: create file paths for hard-coded files
        ############################################################
        finemapping_pips_path = os.path.join(data_dir, "stingseq_data/example_loci/locus_1/finemap_results.txt")
        full_hg19_sumstats = os.path.join(data_dir, "stingseq_data/example_loci/locus_1/variant120-MONO_P-ukbb-variants.hg19.csv")

        alphagenome_gws_scores = os.path.join(data_dir, "stingseq_data/example_loci/locus_1/variant120-MONO_P-ukbb-variants.alphagenome_vep_scores.csv")
        borzoi_prime_gws_scores = os.path.join(data_dir, "stingseq_data/example_loci/locus_1/variant120-MONO_P-ukbb-variants.borzoi_prime_vep_scores.csv")

        hg19_variant_of_interest = "8:56886156" 
        hg19_hit_chrom = hg19_variant_of_interest.split(":")[0]
        hg19_hit_pos = int(hg19_variant_of_interest.split(":")[1])

        ###########################################################
        # step 2: load the supporting data as dictionaries
        ###########################################################

        # create map from hg19 variant to pip
        hg19_variant_to_pip = {}
        pips_df = pd.read_csv(finemapping_pips_path, sep="\s+", dtype=str)

        # note to self: extracting +/- 2 million to be generous 
        pips_df = pips_df[pips_df['chromosome'] == "08"]
        pips_df = pips_df[pips_df['position'].astype(float).sub(hg19_hit_pos).abs() < 2_000_000]

        for _, row in pips_df.iterrows():
            chrom = int(row['chromosome'])
            pos = int(row['position'])
            pip = float(row['prob'])
            assert chrom == int(hg19_hit_chrom) and abs(pos - hg19_hit_pos) <= 2_000_000

            hg19_variant_to_pip[f"{chrom}:{pos}"] = pip
        
        # create map from hg38 variant to alphagenome score
        hg38_variant_to_alphagenome_score = {}
        ag_scores_df = pd.read_csv(alphagenome_gws_scores, dtype=str)

        for _, row in ag_scores_df.iterrows():
            hg38_variant_to_alphagenome_score[f"{row['chrom']}:{row['pos']}"] = float(row['score'])

        # create map from hg38 variant to borzoi prime score
        hg38_variant_to_borzoi_prime_score = {}
        bp_scores_df = pd.read_csv(borzoi_prime_gws_scores, dtype=str)

        for _, row in bp_scores_df.iterrows():
            hg38_variant_to_borzoi_prime_score[f"{row['chrom']}:{row['pos']}"] = float(row['score'])
        
        ###########################################################
        # step 2: iterate through the hg19 sumstats for this locus
        #         and gather all the supplementary data:
        #         - pips from FINEMAP
        #         - alphagenome score (only for gws)
        #         - borzoi prime score (only for gws)
        ###########################################################
        out_fd = open(output.masterfile, "w")
        hg19_sumstats_df = pd.read_csv(full_hg19_sumstats, dtype=str)
        
        out_fd.write("chrom,hg19_pos,hg38_pos,allele1,allele2,pval,-log10p,gws,hit_variant,pip,ag_score,bp_score\n")

        for _, row in hg19_sumstats_df.iterrows():
            chrom = str(row['chrom'])
            pos = int(row['pos'])

            allele1 = row['allele1']
            allele2 = row['allele2']
            pval = row['pval']
            neglog10p = -np.log10(float(pval))

            # define binary description variables
            gws = "Yes" if float(pval) <= 5e-8 else "No"
            hit_status = "Yes" if f"{chrom}:{pos}" == hg19_variant_of_interest else "No"

            # define default values for variables we want to look up
            pip_value = hg19_variant_to_pip.get(f"{chrom}:{pos}", -1)
            ag_score = -1
            bp_score = -1

            # check if we can lift, if not, we will just write it data
            lifted, hg38_chrom, hg38_pos = lift_from_hg19_to_hg38(chrom, pos)
            if not lifted:
                out_fd.write(f"{chrom},{pos},-1,{allele1},{allele2},{pval},{neglog10p},{gws},{hit_status},{pip_value},{ag_score},{bp_score}\n")
                continue 
            
            # we can lift to hg38, so we might have a score ...
            ag_score = hg38_variant_to_alphagenome_score.get(f"{hg38_chrom}:{hg38_pos}", ag_score)
            bp_score = hg38_variant_to_borzoi_prime_score.get(f"{hg38_chrom}:{hg38_pos}", bp_score)

            out_fd.write(f"{chrom},{pos},{hg38_pos},{allele1},{allele2},{pval},{neglog10p},{gws},{hit_status},{pip_value},{ag_score},{bp_score}\n")

