#!/usr/bin/env python3

import argparse
import tensorflow as tf
import tensorflow_hub as hub
import joblib
import gzip
import time
import kipoiseq
from kipoiseq import Interval
import pyfaidx
import pandas as pd
import numpy as np
import enformer_library

ENFORMER_SEQUENCE_CONTEXT = 393216

def one_hot_encode(sequence):
  return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

def compute_l2_diff_log1p_score(ref, alt):
    log_diff = np.log2(1+alt) - np.log2(1+ref)
    l2_norms = np.linalg.norm(log_diff, ord=2, axis=0)    
    return l2_norms

def get_revcomp(seq):
    table = str.maketrans(
        "ACGTacgt",
        "TGCAtgca"
    )
    return seq.translate(table)[::-1]

def run_enformer_scoring(args):
    ################################################
    # section 1: load model and variants
    ################################################
    
    fasta_extractor = enformer_library.FastaStringExtractor(args.genome)
    seq_extractor = kipoiseq.extractors.VariantSeqExtractor(reference_sequence=fasta_extractor)
    
    model = enformer_library.Enformer('https://tfhub.dev/deepmind/enformer/1')
    print(f"\n[log] finished loading the genome and enformer model.")
    
    variants_df = pd.read_csv(args.variant_csv_file, dtype=str)
    print(f"[log] finished loading the variant dataframe.\n")
    
    ######################################################
    # section 2: iterate through variants and score them
    ######################################################
    variant_to_scores = {}
    start_time = time.time()
    
    for index, row in variants_df.iterrows():
        chrom = row['chrom']
        pos = int(row['pos'])
        ref_allele = row['ref']
        alt_allele = row['alt']
        
        print(f"[log] processing chr{chrom}:{pos}:{ref_allele}>{alt_allele} ...")
        
        # extract input context sequence
        variant = kipoiseq.Variant(f"chr{chrom}", pos, ref_allele, alt_allele)
        variant_interval = kipoiseq.Interval(f"chr{chrom}", pos, pos).resize(ENFORMER_SEQUENCE_CONTEXT)
        
        center = variant_interval.center() - variant_interval.start
        ref_seq = seq_extractor.extract(variant_interval, [], anchor=center).upper()
        alt_seq = seq_extractor.extract(variant_interval, [variant], anchor=center).upper()

        # get reverse complement sequences for each
        revcomp_ref_seq = get_revcomp(ref_seq)
        revcomp_alt_seq = get_revcomp(alt_seq)
        
        # make predictions for the reference and alternate allele
        ref_preds = model.predict_on_batch(one_hot_encode(ref_seq)[np.newaxis])['human'][0]
        alt_preds = model.predict_on_batch(one_hot_encode(alt_seq)[np.newaxis])['human'][0]

        revcomp_ref_preds = model.predict_on_batch(one_hot_encode(revcomp_ref_seq)[np.newaxis])['human'][0]
        revcomp_alt_preds = model.predict_on_batch(one_hot_encode(revcomp_alt_seq)[np.newaxis])['human'][0]

        # compute score based on all tracks
        fwd_scores = compute_l2_diff_log1p_score(ref_preds, alt_preds)
        revcomp_scores = compute_l2_diff_log1p_score(revcomp_ref_preds, revcomp_alt_preds)
        
        avg_l2_norms = (fwd_scores+revcomp_scores)/2
        vep_score = np.linalg.norm(avg_l2_norms, ord=2)
        
        print(f"[log] l2 of l2 scores: {vep_score}\n")

        # save score for this variant
        variant_to_scores[f"{chrom}_{pos}_{ref_allele}_{alt_allele}"] = vep_score

    elapsed_time = time.time() - start_time
    sec_per_variant = elapsed_time/variants_df.shape[0]
    
    print(f"[log] finished processing {variants_df.shape[0]} variants ({sec_per_variant:.3f} sec per variant)")
    
    ################################################
    # section 4: write out the scores
    ################################################ 
    with open(args.output, "w") as out_fd:
        out_fd.write("chrom,pos,ref,alt,tool,score\n")
        for var_id, score in variant_to_scores.items():
            chrom, pos, ref, alt = var_id.split("_")
            out_fd.write(f"{chrom},{pos},{ref},{alt},enformer,{score}\n")
    
    print(f"[log] finished writing the variant scores.\n")
        
def parse_arguments():
    parser = argparse.ArgumentParser(description='Run Enformer VEP scoring')
    parser.add_argument('--variant-list', dest="variant_csv_file", required=True, help='Path to CSV file containing variants')
    parser.add_argument('--genome', dest="genome", required=True, help="Path for genome")
    parser.add_argument('--output', dest="output", required=True, help="Path to output file path")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    run_enformer_scoring(args)