######################################################
# Name: analyze_stingseq_gwas_loci.smk
# Description: analyze the vep scores for variants 
#              in sting-seq hit loci
# Date: January 9th, 2026
######################################################

rule analyze_stingseq_hit_loci_percentiles:
    input:
        alphagenome_vep_files=generate_list_of_alphagenome_scores_for_hg38_sumstats(),
        borzoi_vep_files=generate_list_of_borzoi_scores_for_hg38_sumstats(),
        borzoi_prime_vep_files=generate_list_of_borzoi_prime_scores_for_hg38_sumstats(),
        enformer_vep_files=generate_list_of_enformer_scores_for_hg38_sumstats(),
        hit_masterlist="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.hit_loci.csv",
        nonhit_masterlist="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.nonhit_loci.csv"
    output:
        percentiles="analyses/stingseq_datasets/all_gwas_loci/analysis/percentile_scores.csv"
    run:
        import os
        import sys
        import pandas as pd
        import numpy as np

        ##############################################
        # section 0: load in the hit snps and nonhit
        #            snps
        #            - create dictionary mapping
        #              id to snp
        ##############################################        
        hit_id_to_snp = {}
        hit_df = pd.read_csv(input.hit_masterlist, dtype=str)
        for _, row in hit_df.iterrows():
            hit_id_to_snp[int(row['variant_id'])] = f"{row['chrom']}:{row['hg38_pos']}"
        
        nonhit_id_to_snp = {}
        nonhit_df = pd.read_csv(input.nonhit_masterlist, dtype=str)
        for _, row in nonhit_df.iterrows():
            nonhit_id_to_snp[int(row['variant_id'])] = f"{row['chrom']}:{row['hg38_pos']}"
        
        ###############################################
        # section 1: load in each vep score file
        ###############################################
        
        # open output file
        out_fd = open(output.percentiles, "w")
        out_fd.write(f"file,locus_type,tool,chrom,pos,vep_percentile\n")

        # iterate through each locus and compute percentiles
        for tool, vep_file_list in [("alphagenome", input.alphagenome_vep_files), ("enformer", input.enformer_vep_files), ("borzoi", input.borzoi_vep_files), ("borzoi_prime", input.borzoi_prime_vep_files)]: 
            print("#"*20 + f"\n# tool: {tool}\n" + "#"*20)
            for vep_file in vep_file_list:
                print(f"[log] processing: {vep_file}")

                # figure out the locus type and variant id
                loci_type = vep_file.split("/")[-2]
                filename = os.path.basename(vep_file)
                variant_id = int(filename.split("-")[0].replace('variant', ''))

                # determine the lead snp for this locus
                lead_snp = hit_id_to_snp[variant_id] if loci_type == "hit_loci" else nonhit_id_to_snp[variant_id]

                # extract vep scores
                scores_df = pd.read_csv(vep_file, dtype=str)
                all_vep_scores = scores_df['score'].astype(float).tolist()

                # make sure the lead snp is present
                scores_df["snp"] = scores_df["chrom"].astype(str) + ":" + scores_df["pos"].astype(str)
                hit_rows = scores_df[scores_df["snp"].isin([lead_snp])]
                
                if hit_rows.shape[0] != 1: 
                    print(f"[warning] the lead snp for this locus is not present (snp = {lead_snp}, file={vep_file})")
                    continue
                
                # compute the percentile for the lead snp
                hit_score_list = hit_rows['score'].astype(float).tolist()
                hit_percentiles = [np.mean(np.array(all_vep_scores) < val) * 100 for val in hit_score_list]
                
                for index, percentile in enumerate(hit_percentiles):
                    chrom = hit_rows.iloc[index]['chrom']
                    pos = hit_rows.iloc[index]['pos']
                    out_fd.write(f"{filename},{loci_type},{tool},{chrom},{pos},{percentile}\n")


rule analyze_stingseq_loci_size_after_filtering:
    input:
        hit_file="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_counts_per_locus.hit_loci.finalized.csv",
        nonhit_file="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_counts_per_locus.nonhit_loci.finalized.csv"
    output:
        percentiles="analyses/stingseq_datasets/all_gwas_loci/analysis/loci_size.csv"
    run:
        import pandas as pd

        # open output file
        out_fd = open(output.percentiles, "w")
        out_fd.write(f"locus_type,num_variants\n")

        # hard-coded parameters to filter to scored loci
        specific_consort = "ukbb"
        min_variants = 50
        max_num_loci_to_score = 60

        # go through each locus type
        for locus_type, file_path in [('hit_loci', input.hit_file), ('nonhit_loci', input.nonhit_file)]:
            variant_df = pd.read_csv(file_path, dtype=str)

            # subset to relevant sumstats
            variant_df = variant_df[variant_df['consort'].str.contains(specific_consort)]
            variant_df = variant_df[variant_df['num_variants'].astype(int) >= min_variants]
            
            # take the top n loci (to reduce scoring time)
            variant_df = variant_df.head(max_num_loci_to_score)

            for _, row in variant_df.iterrows():
                out_fd.write(f"{locus_type},{row['num_variants']}\n")
        
            
