######################################################
# Name: preprocess_traitgym_datasets.smk
# Description: Preprocesses variant datasets the
#              TraitGym paper
# Date: December 17th, 2025
######################################################

rule process_traitgym_parquet_files:
    output:
        omim_pathogenic="preprocessed_snp_lists/omim_pathogenic_snp_list.csv",
        omim_common="preprocessed_snp_lists/omim_common_snp_list.csv",
        gwas_causal="preprocessed_snp_lists/gwas_causal_snp_list.csv",
        gwas_noncausal="preprocessed_snp_lists/gwas_non_snp_list.csv"
    run:
        import os
        import pandas as pd 

        #############################################
        # step 1: check input files exist 
        #         (present in data_dir)
        #############################################
        gwas_dataset = "snp_datasets/traitgym_gwas_finemapped_nc/test.parquet"
        omim_dataset = "snp_datasets/traitgym_omim_nc/test.parquet"

        for path in [omim_dataset, gwas_dataset]:
            if not os.path.isfile(os.path.join(data_dir, path)):
                print(f"[error] file path is not valid: {os.path.join(data_dir, path)}")

        #############################################
        # step 2: write out the variants 
        #############################################
        def parse_parquet_for_variants(parquet_file, label, outfile):
            df = pd.read_parquet(parquet_file)
            subset_df = df[df['label'] == label]
            
            with open(outfile, "w") as out_fd:
                out_fd.write("chrom,pos,ref,alt\n")
                for index, row in subset_df.iterrows():
                    out_fd.write(f"{row['chrom']},{row['pos']},{row['ref']},{row['alt']}\n")
        
        parse_parquet_for_variants(os.path.join(data_dir, omim_dataset),
                                   True,
                                   output.omim_pathogenic)
        parse_parquet_for_variants(os.path.join(data_dir, omim_dataset),
                                   False,
                                   output.omim_common)
        parse_parquet_for_variants(os.path.join(data_dir, gwas_dataset),
                                   True,
                                   output.gwas_causal)
        parse_parquet_for_variants(os.path.join(data_dir, gwas_dataset),
                                   False,
                                   output.gwas_noncausal)
        