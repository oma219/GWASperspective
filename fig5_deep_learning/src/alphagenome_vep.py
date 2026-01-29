
import os
import time
import argparse
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

def score_variants(args):
    ##############################################
    # step 1: load in the variants
    ##############################################
    variants_df = pd.read_csv(args.variant_list)    
    print(f"\n[log] loaded dataframe with n = {variants_df.shape[0]} variants.\n")

    ###############################################################
    # step 2: initialize alphagenome model for querying with API
    ###############################################################
    dna_model = dna_client.create(args.api_key)    
    sequence_length = args.seq_length
    
    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[f"SEQUENCE_LENGTH_{sequence_length}"]
    organism = dna_client.Organism.HOMO_SAPIENS
    
    ######################################
    # Note: choose these tracks to match
    #       TraitGym's benchmarking
    ######################################
    tracks = [
        "ATAC",
        "DNASE",
        "CHIP_TF",
        "CHIP_HISTONE",
        "CAGE",
        "PROCAP",
        "RNA_SEQ",
    ]
    scorers = [
        variant_scorers.CenterMaskScorer(
            requested_output=getattr(dna_client.OutputType, track),
            width=None,
            aggregation_type=variant_scorers.AggregationType.L2_DIFF_LOG1P,
        )
        for track in tracks
    ]

    #######################################
    # step 3: run scoring for all variants
    #######################################
    variant_to_scores = {}
    start_time = time.time()
    
    for index, row in variants_df.iterrows():
        chrom = f"chr{row['chrom']}"
        pos = int(row['pos'])
        ref = str(row['ref'])
        alt = str(row['alt'])
        print(f"[log] processing variant = {chrom}_{pos}_{ref}_{alt} ... ")

        try:
            # create variant and interval object
            variant = genome.Variant(
                chromosome=chrom,
                position=pos,
                reference_bases=ref,
                alternate_bases=alt
            )
            interval = variant.reference_interval.resize(sequence_length)
            
            # create a forward/reverse complement strand
            interval_fwd = interval.copy()
            interval_fwd.strand = "+"
            interval_revcomp = interval.copy()
            interval_revcomp.strand = "-"
            
            def score_given_interval(interval):
                variant_scores = dna_model.score_variant(
                        interval=interval,
                        variant=variant,
                        variant_scorers=scorers,
                        organism=organism
                )
                return variant_scorers.tidy_scores([variant_scores])
            
            # compute scores for both strands
            fwd_scores = score_given_interval(interval_fwd)
            revcomp_scores = score_given_interval(interval_revcomp)
            
            # make sure they have the same columns
            assert (fwd_scores.index == revcomp_scores.index).all() and (fwd_scores.columns == revcomp_scores.columns).all()

            # compute average score across strands
            scores_df = fwd_scores.copy()
            scores_df['avg_raw_score'] = (fwd_scores.raw_score + revcomp_scores.raw_score)/2
            
            max_l2_norm_score = scores_df['avg_raw_score'].max()
            print(f"[log] max of l2 norms: {max_l2_norm_score}\n")
            
            # save score for this variant
            variant_to_scores[f"{chrom.replace('chr', '')}_{pos}_{ref}_{alt}"] = max_l2_norm_score
            
        except Exception as e:
            print(f"[error] {e}"); exit(1)
    
    elapsed_time = time.time() - start_time
    sec_per_variant = elapsed_time/variants_df.shape[0]
    
    print(f"[log] finished processing {variants_df.shape[0]} variants ({sec_per_variant:.3f} sec per variant)")
    
    #######################################
    # step 4: write out the scores
    #######################################
    with open(args.output, "w") as out_fd:
        out_fd.write("chrom,pos,ref,alt,tool,score\n")
        for var_id, score in variant_to_scores.items():
            chrom, pos, ref, alt = var_id.split("_")
            out_fd.write(f"{chrom},{pos},{ref},{alt},alphagenome,{score}\n")
    
    print(f"[log] finished writing the variant scores.\n")
    

def parse_arguments():
    parser = argparse.ArgumentParser(description="Score variants using AlphaGenome and save median assay values per variant")
    parser.add_argument('--variant-list', dest="variant_list", required=True, help='Path to the variant list CSV file')
    parser.add_argument('--api-key', dest="api_key", required=True, help="API key for accessing AlphaGenome API")
    parser.add_argument('--seq-length', dest="seq_length", required=True, help="Short string for size of AlphaGenome context (e.g. 1MB for 1048576)")
    parser.add_argument('--output', dest="output", required=True, help="Path to output file path")
    
    args = parser.parse_args()
    
    ############################################
    # check 1: check variant list looks valid
    ############################################
    if not os.path.exists(args.variant_list):
        print("\n[error] variant list path is not valid.\n"); exit(1)
        
    variants_df = pd.read_csv(args.variant_list)
    if not set(["chrom", "pos", "ref", "alt"]).issubset(variants_df.columns):
        print("[error] csv must contain chrom, pos, ref, alt columns\n"); exit(1)
    
    ###########################################
    # check 2: make sure sequence length is 
    #          a valid option
    ###########################################
    if args.seq_length not in ["2KB", "16KB", "100KB", "500KB", "1MB"]:
        print(f"[error] provided sequence length is not valid: {args.seq_length}\n");
    
    return args

if __name__ == "__main__":
    args = parse_arguments()
    score_variants(args)