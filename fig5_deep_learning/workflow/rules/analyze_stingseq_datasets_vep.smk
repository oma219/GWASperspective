######################################################
# Name: analyze_traitgym_datasets_vep.smk
# Description: analyze the vep scores for variants
#              from the traitgym datasets
# Date: January 3rd, 2026
######################################################

rule analyze_stingseq_dataset_vep:
    input:
        vep_files=expand("variant_scores/{tool}/{dataset}_{tool}_vep_scores.csv",
                         tool=['alphagenome', 'borzoi', 'borzoi_prime', 'enformer', 'chrombpnet', 'gpn_msa'],
                         dataset=['stingseq_hits', 'stingseq_nonhits'])
    output:
        analysis_csv="analyses/stingseq_datasets/binary_classification/auprc_on_stingseq_datasets.csv",
        pr_curve_csv="analyses/stingseq_datasets/binary_classification/precision_recall_curve_on_stingseq_datasets.csv"
    run:
        import os
        import pandas as pd
        from sklearn.metrics import average_precision_score
        from sklearn.metrics import precision_recall_curve

        ##############################################
        # section 0: define helper methods
        ##############################################
        def compute_auprc_from_df(input_df, pos_class, neg_class, tool):
            pos_class_scores = input_df[(input_df['dataset'] == pos_class) & (input_df['tool'] == tool)]['score'].tolist()
            neg_class_scores = input_df[(input_df['dataset'] == neg_class) & (input_df['tool'] == tool)]['score'].tolist()

            pos_class_scores = [float(x) for x in pos_class_scores]
            neg_class_scores = [float(x) for x in neg_class_scores]

            all_scores = pos_class_scores + neg_class_scores
            all_labels = [1 for x in range(len(pos_class_scores))] + [0 for x in range(len(neg_class_scores))]
            auprc = round(average_precision_score(all_labels, all_scores), 4)
            precision, recall, thresholds = precision_recall_curve(all_labels, all_scores)
            return auprc, precision, recall

        ##############################################
        # section 1: read all individual vep files
        #            combine into 1 file
        ##############################################
        all_dfs = []
        for curr_file in input.vep_files:
            curr_df = pd.read_csv(curr_file, dtype=str)
            curr_df["score"] = pd.to_numeric(curr_df["score"], errors="raise")
            assert list(curr_df.columns[:6]) == ['chrom', 'pos', 'ref', 'alt', 'tool', 'score'], f"Columns are: {curr_df.columns.tolist()}"

            file_basename = os.path.basename(curr_file)
            curr_dataset = "_".join(file_basename.split("_")[:2])
            
            curr_df.insert(5, 'dataset', curr_dataset)

            # take maximum score for snps where we score multiple alt alleles
            compressed_curr_df = curr_df.loc[curr_df.groupby(["chrom", "pos"])["score"].idxmax()].reset_index(drop=True)
            
            all_dfs.append(compressed_curr_df)
        
        combined_df = pd.concat(all_dfs, axis=0, ignore_index=True)

        print(f"\n[log] total number of rows after combining: {combined_df.shape[0]}")
        print(f"[log] unique tools: {combined_df['tool'].unique()}\n")

        ##############################################
        # section 2: open output file for auprcs
        ##############################################
        auprc_fd = open(output.analysis_csv, "w")
        auprc_fd.write("dataset,tool,author,auprc\n")

        pr_fd = open(output.pr_curve_csv, "w")
        pr_fd.write("dataset,tool,precision,recall\n")

        ##############################################
        # section 3: iterate through all the tools
        #            and datasets combos and compute
        #            in-house auprcs
        ##############################################
        for dataset, pos_class, neg_class in [('stingseq', 'stingseq_hits', 'stingseq_nonhits')]:
            for tool in combined_df['tool'].unique():
                auprc, precision, recall = compute_auprc_from_df(combined_df, pos_class, neg_class, tool)
                auprc_fd.write(f"{dataset},{tool},homebrew,{auprc}\n")

                for x, y in zip(precision, recall):
                    pr_fd.write(f"{dataset},{tool},{x},{y}\n")

rule analyze_stingseq_dataset_vep_across_dist_and_functional_marks:
    input:
        vep_files=expand("variant_scores/{tool}/{dataset}_{tool}_vep_scores.csv",
                         tool=['alphagenome', 'borzoi', 'borzoi_prime', 'enformer', 'chrombpnet', 'gpn_msa'],
                         dataset=['stingseq_hits', 'stingseq_nonhits']),
        hits_list="preprocessed_snp_lists/stingseq_hits_snp_list.csv",
        nonhits_list="preprocessed_snp_lists/stingseq_nonhits_snp_list.csv"
    output:
        analysis_csv="analyses/stingseq_datasets/binary_classification/stingseq_variant_scores_with_annotations.csv"
    run:
        import os
        import pandas as pd
        import numpy as np
        from liftover import get_lifter
        
        converter = get_lifter('hg19', 'hg38', one_based=True)

        ##################################################
        # section 0: define helper methods
        ##################################################
        def lift_from_hg19_to_hg38(hg19_chrom, hg19_pos):
            pos_lifted = converter[f"chr{hg19_chrom}"][int(hg19_pos)]
            if len(pos_lifted) == 0:
                return False, "", ""
            else:
                hg38_chrom = pos_lifted[0][0].replace("chr", "")
                hg38_pos = int(pos_lifted[0][1])
                assert hg38_chrom == hg19_chrom

                return True, hg38_chrom, hg38_pos

        def robust_minmax(x, lower=0.01, upper=0.99):
            lo = x.quantile(lower)
            hi = x.quantile(upper)
            return ((x - lo) / (hi - lo)).clip(0, 1)

        ######################################################
        # section 1: load in a hit/non-hit sets of variants
        ######################################################
        hit_set = set(); nonhit_set = set()

        hit_df = pd.read_csv(input.hits_list, dtype=str)
        for index, row in hit_df.iterrows():
            hit_set.add(f"{row['chrom']}:{row['pos']}")

        nonhit_df = pd.read_csv(input.nonhits_list, dtype=str)
        for index, row in nonhit_df.iterrows():
            nonhit_set.add(f"{row['chrom']}:{row['pos']}")
        
        print(f"\n[log] number of hit variants: {len(hit_set)}")
        print(f"[log] number of non-hit variants: {len(nonhit_set)}\n")

        ######################################################
        # section 2: load in variant annotations 
        #            - absolulte dist to TSS of correct gene 
        #              (closest one if multiple, -1 if there
        #               is no target gene)
        #            - functional marks score 
        #######################################################
        snp_to_annot = {}

        full_variant_df = pd.read_csv(os.path.join(data_dir, "stingseq_data/variant_list/john_table3f_v3.csv"), dtype=str)
        for index, row in full_variant_df.iterrows():
            snp = row['SNP Coordinates (hg19)']
            hg19_chrom = snp.split(":")[0]
            hg19_pos = snp.split(":")[1]

            lifted, hg38_chrom, hg38_pos = lift_from_hg19_to_hg38(hg19_chrom, hg19_pos)

            # move on if this variant cannot be lifted
            if not lifted:
                continue

            snp = f"{hg38_chrom}:{hg38_pos}"
            if snp not in snp_to_annot:
                snp_to_annot[snp] = {"tss_dist": -1, "func_score": int(row['Functional Score'])}

            # check if this row corresponds to a target gene
            if pd.notna(row["Q-value (500 kb)"]) and (row["Q-value (500 kb)"] != "") and (float(row["Q-value (500 kb)"]) < 0.05) and snp_to_annot[snp]["tss_dist"] != -1:
                snp_to_annot[snp]["tss_dist"] = min(abs(int(row['TSS Distance'])), snp_to_annot[snp]["tss_dist"])
            elif pd.notna(row["Q-value (500 kb)"]) and (row["Q-value (500 kb)"] != "") and (float(row["Q-value (500 kb)"]) < 0.05) and snp_to_annot[snp]["tss_dist"] == -1:
                snp_to_annot[snp]["tss_dist"] = abs(int(row['TSS Distance']))

        ######################################################
        # section 3: load in all VEP scores for all variants 
        ######################################################
        snp_to_vep_scores = {}

        for tool in ['alphagenome', 'borzoi', 'borzoi_prime', 'enformer', 'chrombpnet', 'gpn_msa']:
            for dataset in ['stingseq_hits', 'stingseq_nonhits']:

                vep_file_path = f"variant_scores/{tool}/{dataset}_{tool}_vep_scores.csv"
                vep_df = pd.read_csv(vep_file_path, dtype=str)

                for index, row in vep_df.iterrows():
                    snp = f"{row['chrom']}:{row['pos']}"
                    if snp not in snp_to_vep_scores:
                        snp_to_vep_scores[snp] = {}
                    snp_to_vep_scores[snp][tool] = float(row['score'])
        
        ######################################################
        # section 4: print out the scores including all
        #            functional annotations
        ######################################################
        out_fd = open(output.analysis_csv, "w")
        out_fd.write("chrom,pos,class,dist_to_tss,dist_group,func_score,func_group,tool,score\n")
        
        for variant_set, set_name in [(hit_set, "stingseq_hits"), (nonhit_set, "stingseq_nonhits")]:
            
            dist_group_labels = []; func_group_labels = []
            print(f"[log] processing {set_name} variants (n = {len(variant_set)}) ...")

            for snp in variant_set:
                chrom = snp.split(":")[0]
                pos = snp.split(":")[1]

                dist_to_tss = snp_to_annot[snp]['tss_dist']
                func_score = snp_to_annot[snp]['func_score']

                # assign the variant to distance group
                if dist_to_tss < 10_000:
                    dist_group = 1
                elif 10_000 <= dist_to_tss < 50_000:
                    dist_group = 2
                else:
                    dist_group = 3
                dist_group_labels.append(dist_group)

                # assign variant to functional group
                if func_score == 7:
                    func_group = 1
                elif func_score in [6, 5, 3]:
                    func_group = 2
                else:
                    func_group = 3
                func_group_labels.append(func_group)

                # write out the scores for this variant for each tool
                for tool, score in snp_to_vep_scores[snp].items():
                    out_fd.write(f"{chrom},{pos},{set_name},{dist_to_tss},{dist_group},{func_score},{func_group},{tool},{score}\n")

            dist_group_counts = {label: dist_group_labels.count(label) for label in set(dist_group_labels)}
            func_group_counts = {label: func_group_labels.count(label) for label in set(func_group_labels)}

            print(f"[log] breakdown of different distance groupings: {dist_group_counts}")
            print(f"[log] breakdown of different functional groupings: {func_group_counts}\n")

        print(f"[log] finished writing the full score file")
        out_fd.close()

        ######################################################
        # section 5: load the full score file and normalize
        #            all scores within-tool
        ######################################################
        score_df = pd.read_csv(output.analysis_csv, dtype=str)
        score_df["normalized_score"] = (score_df.astype({"score": float}).groupby("tool", group_keys=False)["score"].apply(robust_minmax))

        score_df.to_csv(output.analysis_csv, index=False)
        print("[log] computed normalized scores and saved to file\n")