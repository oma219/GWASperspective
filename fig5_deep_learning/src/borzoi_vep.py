#!/usr/bin/env python

from __future__ import print_function

from optparse import OptionParser
import json
import pdb
import pickle
import os
import sys
import time
import argparse

import h5py
import numpy as np
import pandas as pd
import pysam
from scipy.sparse import dok_matrix
from scipy.special import rel_entr
import tensorflow as tf
from tqdm import tqdm

import kipoiseq
from kipoiseq import Interval
import pyfaidx

from baskerville import dataset
from baskerville import seqnn
from baskerville import stream
from baskerville import vcf as bvcf

class FastaStringExtractor:
    
    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                          trimmed_interval.start + 1,
                                          trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()

def one_hot_encode(sequence):
  return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

def get_revcomp(seq):
    table = str.maketrans(
        "ACGTacgt",
        "TGCAtgca"
    )
    return seq.translate(table)[::-1]

def compute_l2_diff_log1p_scores(ref, alt):
    log_diff = np.log2(1+ref) - np.log2(1+alt)
    l2_norms = np.linalg.norm(log_diff, axis=0)    
    return l2_norms

def score_variants(args):
    #################################################
    # section 1: load model parameters and targets 
    #################################################
    with open(args.model_params) as params_fd:
        params = json.load(params_fd)
        
    params_model = params['model']
    params_train = params['train']
    
    ################################################
    # section 2: load model and genome
    ################################################
    seqnn_model = seqnn.SeqNN(params_model)
    seqnn_model.restore(args.model)
    
    seqnn_model.build_slice(None) # To-do: optionally specific specific tracks to output

    input_seq_length = params_model['seq_length']
    targets_length = seqnn_model.target_lengths[0]
    num_targets = seqnn_model.num_targets()
    
    print(f"\n[log] size of input context: {input_seq_length}")
    print(f"[log] number of genomic bins in each track output (32 bp each): {targets_length}")
    print(f"[log] number of output tracks: {num_targets}\n")
    
    genome = pysam.Fastafile(args.genome)
    print(f"[log] loaded this genome for extracting context: {os.path.basename(args.genome)}\n")
    
    fasta_extractor = FastaStringExtractor(args.genome)
    seq_extractor = kipoiseq.extractors.VariantSeqExtractor(reference_sequence=fasta_extractor)
    
    ################################################
    # section 3: iterate through the variants
    ################################################ 
    variant_to_scores = {}
    start_time = time.time()
    
    variants_df = pd.read_csv(args.variant_list, dtype=str)
    for index, row in variants_df.iterrows():
        chrom = row['chrom']
        pos = int(row['pos'])
        ref = row['ref']
        alt = row['alt']
        
        print(f"[log] processing variant = {chrom}:{pos}:{ref}>{alt} ...")

        # extract the context sequences
        variant = kipoiseq.Variant(f"chr{chrom}", pos, ref, alt)
        variant_interval = kipoiseq.Interval(f"chr{chrom}", pos, pos).resize(input_seq_length)
        
        center = variant_interval.center() - variant_interval.start
        ref_seq = seq_extractor.extract(variant_interval, [], anchor=center)
        alt_seq = seq_extractor.extract(variant_interval, [variant], anchor=center)
        
        # get reverse complement sequences for each
        revcomp_ref_seq = get_revcomp(ref_seq)
        revcomp_alt_seq = get_revcomp(alt_seq)
                
        # make predictions for the reference and alternate allele (as well as revcomp)
        ref_preds = seqnn_model(one_hot_encode(ref_seq)[np.newaxis])[0]
        alt_preds = seqnn_model(one_hot_encode(alt_seq)[np.newaxis])[0]     
        
        revcomp_ref_seq = seqnn_model(one_hot_encode(revcomp_ref_seq)[np.newaxis])[0]
        revcomp_alt_seq = seqnn_model(one_hot_encode(revcomp_alt_seq)[np.newaxis])[0]
        
        # summarize the output tracks
        fwd_scores = compute_l2_diff_log1p_scores(ref_preds, alt_preds)
        revcomp_scores = compute_l2_diff_log1p_scores(revcomp_ref_seq, revcomp_alt_seq)
        
        avg_l2_norms = (fwd_scores+revcomp_scores)/2
        vep_score = np.linalg.norm(avg_l2_norms)
        
        print(f"[log] l2 of l2 scores: {vep_score}\n")

        # save score for this variant
        variant_to_scores[f"{chrom}_{pos}_{ref}_{alt}"] = vep_score

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
            out_fd.write(f"{chrom},{pos},{ref},{alt},{args.tool_name},{score}\n")
    
    print(f"[log] finished writing the variant scores.\n")
 
def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to run Borzoi model for VEP")
    parser.add_argument('--model-params', dest="model_params", type=str, required=True, help='Path to the model parameters JSON file')
    parser.add_argument('--model', dest="model", type=str, required=True, help='Path to the trained model weights (H5 file)')
    parser.add_argument('--variant-list', dest="variant_list", type=str, required=True, help='CSV file containing variants to score')
    parser.add_argument('--genome', dest="genome", type=str, required=True, help='Path to the genome FASTA file')
    parser.add_argument('--output', dest="output", required=True, help="Path to output file path")
    parser.add_argument('--tool-name', default="borzoi", dest="tool_name", help="Name of tool to be used in the output file")
    args = parser.parse_args()
    
    # check #1: make sure all file paths exist
    for path in [args.model_params, args.model, args.variant_list]:
        if not os.path.isfile(path):
            print(f"[error] file path is not valid: {path}"); exit(1)
    
    # check 2: make sure the variant dataframe contains the needed columns
    variants_df = pd.read_csv(args.variant_list, dtype=str)
    if not set(["chrom", "pos", "ref", "alt"]).issubset(variants_df.columns):
        print(f"\n[error] missing columns in variant list\n"); exit(1)
    
    return args

if __name__ == "__main__":
    args = parse_arguments()
    score_variants(args)
    