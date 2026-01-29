######################################################
# Name: analyze_traitgym_datasets_vep.smk
# Description: analyze the vep scores for variants
#              from the traitgym datasets
# Date: January 3rd, 2026
######################################################

rule analyze_traitgym_dataset_vep:
    input:
        vep_files=expand("variant_scores/{tool}/{dataset}_{tool}_vep_scores.csv",
                         tool=['alphagenome', 'borzoi', 'borzoi_prime', 'enformer', 'chrombpnet', 'gpn_msa'],
                         dataset=['omim_pathogenic', 'omim_common', 'gwas_causal', 'gwas_non'])
    output:
        analysis_csv="analyses/traitgym_datasets/auprc_on_traitgym_datasets.csv"
    run:
        import os
        import pandas as pd
        from sklearn.metrics import average_precision_score

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
            return auprc
        
        def compute_auprc_from_traitgym_predictions(dataset, tool):
            variant_labels = []; variant_scores = []

            variant_df = pd.read_parquet(os.path.join(data_dir, f"traitgym_data/predictions/{dataset}/variant_list.parquet"))
            score_df = pd.read_parquet(os.path.join(data_dir, f"traitgym_data/predictions/{dataset}/{tool}.parquet"))

            variant_labels = [bool(label) for label in variant_df['label']]
            variant_scores = [abs(float(score)) for score in score_df['score']]

            group1_scores = [score for label, score in zip(variant_labels, variant_scores) if label]
            group2_scores = [score for label, score in zip(variant_labels, variant_scores) if not label]

            all_scores = group1_scores + group2_scores
            all_labels = [1 for x in range(len(group1_scores))] + [0 for x in range(len(group2_scores))]
            auprc = round(average_precision_score(all_labels, all_scores), 4)
            return auprc

        ##############################################
        # section 1: read all individual vep files
        #            combine into 1 file
        ##############################################
        all_dfs = []
        for curr_file in input.vep_files:
            curr_df = pd.read_csv(curr_file, dtype=str)
            assert list(curr_df.columns[:6]) == ['chrom', 'pos', 'ref', 'alt', 'tool', 'score'], f"Columns are: {curr_df.columns.tolist()}"

            file_basename = os.path.basename(curr_file)
            curr_dataset = "_".join(file_basename.split("_")[:2])
            
            curr_df.insert(5, 'dataset', curr_dataset)
            all_dfs.append(curr_df)
        
        combined_df = pd.concat(all_dfs, axis=0, ignore_index=True)

        print(f"\n[log] total number of rows after combining: {combined_df.shape[0]}")
        print(f"[log] unique tools: {combined_df['tool'].unique()}\n")

        ##############################################
        # section 2: open output file for auprcs
        ##############################################
        out_fd = open(output.analysis_csv, "w")
        out_fd.write("dataset,tool,author,auprc\n")

        ##############################################
        # section 3: iterate through all the tools
        #            and datasets combos and compute
        #            in-house auprcs
        ##############################################
        for dataset, pos_class, neg_class in [('omim', 'omim_pathogenic', 'omim_common'), ('gwas', 'gwas_causal', 'gwas_non')]:
            for tool in combined_df['tool'].unique():
                auprc = compute_auprc_from_df(combined_df, pos_class, neg_class, tool)
                out_fd.write(f"{dataset},{tool},homebrew,{auprc}\n")
        
        ########################################
        # section 4: compute auprcs using the 
        #            predictions from traitgym
        ########################################
        for dataset in ['omim', 'gwas']:
            for tool in ['alphagenome', 'borzoi', 'enformer', 'gpn_msa']:
                auprc = compute_auprc_from_traitgym_predictions(dataset, tool)
                out_fd.write(f"{dataset},{tool},traitgym,{auprc}\n")
        
        out_fd.close()


