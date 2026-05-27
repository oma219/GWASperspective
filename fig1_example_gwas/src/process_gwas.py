#!/usr/bin/env python3

#############################################
# ./process_gwas.py --input-folder /path/
#############################################

import os
import sys
import argparse 
import pandas as pd

def main(args):
    ############################################
    # step 1: identify all unique loci
    ############################################
    leadsnp_to_fm_output = {}
    susie_files = [os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir) if f.endswith('.tsv')]
    
    for curr_file in susie_files:
        leadsnp = os.path.basename(curr_file).split("_")[1]
        if leadsnp not in leadsnp_to_fm_output:
            leadsnp_to_fm_output[leadsnp] = curr_file
    
    print(f"\n[log] number of loci fine-mapped for this gwas: {len(leadsnp_to_fm_output)}\n")
    
    ############################################
    # step 2: compute statistics
    #         - number of loci with at least
    #           1 variant with PIP >= 0.10
    #         - number of variants between
    #           0.01-1.00
    #         - number of variants above 
    #           different thresholds
    ############################################
    num_loci_with_large_pip = 0; num_loci_with_no_large_pips = 0
    num_variants_above_threshold = {0.01: 0,
                                    0.1: 0,
                                    0.5: 0,
                                    0.8: 0,
                                    0.9: 0,
                                    0.99: 0}
    num_variants_in_range = {"0.01-0.1": 0,
                             "0.1-0.5": 0,
                             "0.5-0.8": 0,
                             "0.8-0.9": 0,
                             "0.9-0.99": 0,
                             "0.99-1.0": 0}

    for lead_snp, fm_file in leadsnp_to_fm_output.items():
        df = pd.read_csv(fm_file, sep='\t', dtype=str)
        
        # first, keep track of how many variants above certain thresholds
        num_variants_above_threshold1 = 0
        num_variants_above_threshold2 = 0
        num_variants_above_threshold3 = 0
        num_variants_above_threshold4 = 0
        num_variants_above_threshold5 = 0
        num_variants_above_threshold6 = 0
        
        num_variants_above_threshold1 = df[df['PIP'].astype(float) >= 0.01].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        num_variants_above_threshold2 = df[df['PIP'].astype(float) >= 0.1].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        num_variants_above_threshold3 = df[df['PIP'].astype(float) >= 0.50].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        num_variants_above_threshold4 = df[df['PIP'].astype(float) >= 0.80].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        num_variants_above_threshold5 = df[df['PIP'].astype(float) >= 0.90].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        num_variants_above_threshold6 = df[df['PIP'].astype(float) >= 0.99].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        
        num_variants_above_threshold[0.01] += num_variants_above_threshold1 
        num_variants_above_threshold[0.1] += num_variants_above_threshold2 
        num_variants_above_threshold[0.5] += num_variants_above_threshold3 
        num_variants_above_threshold[0.8] += num_variants_above_threshold4 
        num_variants_above_threshold[0.9] += num_variants_above_threshold5 
        num_variants_above_threshold[0.99] += num_variants_above_threshold6 
        
        # second, document if locus contains a 'high' pip
        if num_variants_above_threshold2 > 0:
            num_loci_with_large_pip += 1
        else:
            num_loci_with_no_large_pips += 1
        
        # third, document number of variants in different ranges
        num_variants_in_range1 = 0
        num_variants_in_range2 = 0
        num_variants_in_range3 = 0
        num_variants_in_range4 = 0
        num_variants_in_range5 = 0
        num_variants_in_range6 = 0

        num_variants_in_range1 = df[(df['PIP'].astype(float) >= 0.01) & (df['PIP'].astype(float) < 0.1)].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        num_variants_in_range2 = df[(df['PIP'].astype(float) >= 0.1) & (df['PIP'].astype(float) < 0.5)].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        num_variants_in_range3 = df[(df['PIP'].astype(float) >= 0.5) & (df['PIP'].astype(float) < 0.8)].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        num_variants_in_range4 = df[(df['PIP'].astype(float) >= 0.8) & (df['PIP'].astype(float) < 0.9)].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        num_variants_in_range5 = df[(df['PIP'].astype(float) >= 0.9) & (df['PIP'].astype(float) < 0.99)].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        num_variants_in_range6 = df[(df['PIP'].astype(float) >= 0.99)].drop_duplicates(subset=['chromosome', 'position']).shape[0]
        
        num_variants_in_range["0.01-0.1"] += num_variants_in_range1
        num_variants_in_range["0.1-0.5"] += num_variants_in_range2
        num_variants_in_range["0.5-0.8"] += num_variants_in_range3
        num_variants_in_range["0.8-0.9"] += num_variants_in_range4
        num_variants_in_range["0.9-0.99"] += num_variants_in_range5
        num_variants_in_range["0.99-1.0"] += num_variants_in_range6

    ############################################
    # step 3: write out statistics
    ############################################
    print(f"[log] number of loci with at lesat 1 variant with PIP >= 0.10: {num_loci_with_large_pip}")
    print(f"[log] number of loci without at least 1 variant with PIP >= 0.10: {num_loci_with_no_large_pips}\n")
    
    print(f"[log] number of variants in different ranges:\n")
    
    print(f"low,high,num_variants")
    for pip_range in ["0.01-0.1", "0.1-0.5", "0.5-0.8", "0.8-0.9", "0.9-0.99", "0.99-1.0"]:
        low = pip_range.split("-")[0]
        high = pip_range.split("-")[1]
        num_variants = num_variants_in_range[pip_range]
        print(f"{low},{high},{num_variants}")

    print(f"\n[log] number of variants above different thresholds:\n")
    
    print(f"threshold,num_variants,num_crispr_guides,num_mpra_alleles,num_cells_for_moi_1_fifth,num_cells_for_moi_1, num_cells_for_moi_5")
    for threshold in [0.01, 0.1, 0.5, 0.8, 0.9, 0.99]:
        num_variants = num_variants_above_threshold[threshold]
        print(f"{threshold},{num_variants},{num_variants * 4},{num_variants * 50},{num_variants * 4 * 100/0.2},{num_variants * 4 * 100/1},{num_variants * 4 * 100/5}")
         

    
def parse_arguments():
    parser = argparse.ArgumentParser(description="Process GWAS input files.")
    parser.add_argument('--input-folder', dest="input_dir", type=str, required=True, help='Path to the input folder')
    args = parser.parse_args()
    
    if not os.path.isdir(args.input_dir):
        print(f"\n[error] input folder is not a valid directory: {args.input_dir}")
        exit(1)
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    main(args)