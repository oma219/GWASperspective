######################################################
# Name: preprocess_stingseq_gwas_loci.smk
# Description: preprocess all the loci from the 
#              sting-seq screen where a hit was 
#              identified to prepare them for 
#              scoring for deep learning tools.
#
#              the difficulty comes from the fact 
#              that different hits come from different
#              traits therefore, it takes some 
#              detective work to trace it back to
#              the correct sumstats file
# Date: january 6th, 2025
######################################################

#######################################################
# section 0: preprocess the data to prepare for 
#            extracting variants from hit loci and
#            nonhit loci
#######################################################

rule identify_traits_in_each_gwas_set:
    output:
        trait_list="analyses/stingseq_datasets/all_gwas_loci/metadata/trait_to_consortium.csv"
    run:
        import os
        import sys

        ukbb_sumstats_dir = os.path.join(data_dir, "stingseq_data/gwas/ukbb")
        bcc_sumstats_dir = os.path.join(data_dir, "stingseq_data/gwas/bcc")

        ukbb_traits = [file.replace('.tsv', '') for file in os.listdir(ukbb_sumstats_dir) if file.endswith('.tsv')]
        bcc_traits = [file.replace('.tsv.gz', '') for file in os.listdir(bcc_sumstats_dir) if file.endswith('.tsv.gz')]

        all_traits_set = set(ukbb_traits + bcc_traits)

        with open(output.trait_list, "w") as out_fd:
            out_fd.write(f"trait,consortium_list\n")
            for trait in all_traits_set:
                if trait in ukbb_traits and trait in bcc_traits:
                    out_fd.write(f"{trait},ukbb|bcc\n")
                elif trait in ukbb_traits:
                    out_fd.write(f"{trait},ukbb\n")
                elif trait in bcc_traits:
                    out_fd.write(f"{trait},bcc\n")
                else:
                    print(f"[error] trait should be present in either ukbb or bcc lists.")
                    exit(1)
        
        print(f"\n[log] finished writing the trait-to-consoritum list.\n")


rule build_variant_masterlist_for_all_hit_loci_analysis:
    input:
        hit_list="preprocessed_snp_lists/stingseq_hits_snp_list.csv",
        trait_list="analyses/stingseq_datasets/all_gwas_loci/metadata/trait_to_consortium.csv"
    output:
        masterlist="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.hit_loci.csv"
    run:
        import os
        import sys
        import pandas as pd 
        from liftover import get_lifter
        
        converter = get_lifter('hg19', 'hg38', one_based=True)
        ###########################################################
        # section 0: load in various dictionaries
        #            - trait-to-consortiums
        #            - hg19_variant-to-trait
        ###########################################################

        # load in a list of traits and what consortiums has those gwas
        trait_df = pd.read_csv(input.trait_list)
        trait_to_consortium = {}

        for _, row in trait_df.iterrows():
            trait_to_consortium[row['trait']] = row['consortium_list'].split("|")
        
        print(f"\n[log] loaded in list of traits (n = {len(trait_to_consortium)} unique traits)")

        # load in variant metadata to map variant to its gwas trait
        variants_df = pd.read_csv(os.path.join(data_dir, "stingseq_data/variant_list/john_table3f_v3.csv"), dtype=str)
        hg19_variant_to_trait = {}

        # cases where the variant file has a gwas trait that doesn't match table S1A
        special_cases = {"EO": "EOS",
                         "LYMPH_P": "LYM_P",
                         "NEUT": "NEU",
                         "EO_P": "EOS_P",
                         "LYMPH": "LYM",
                         "NEUT_P": "NEU_P",
                         "MONO": "MON",
                         "RBC_W": "RDW", # verified this one by hand
                         "BASO_P": "BAS_P"}

        # creates the variant to trait map
        for _, row in variants_df.iterrows():
            top_trait = row['Top GWAS fine-mapped']
            assert top_trait in special_cases or top_trait in trait_to_consortium, f"[error] unexpected trait found = {top_trait}"
            if top_trait in special_cases:
                top_trait = special_cases[top_trait]
            
            hg19_variant_to_trait[row['SNP Coordinates (hg19)']] = top_trait
            
        print("[log] finished building the dictionaries\n")

        ###########################################################
        # section 1: write variant masterlist
        ###########################################################

        def lift_from_hg38_to_hg19(hg38_chrom, hg38_pos):
            converter = get_lifter('hg38', 'hg19', one_based=True)
            pos_lifted = converter[f"chr{hg38_chrom}"][int(hg38_pos)]

            if len(pos_lifted) == 0:
                return False, "", ""
            else:
                hg19_chrom = pos_lifted[0][0].replace("chr", "")
                hg19_pos = int(pos_lifted[0][1])
                assert hg19_chrom == hg38_chrom

                return True, hg19_chrom, hg19_pos

        with open(output.masterlist, "w") as out_fd:
            out_fd.write("variant_id,chrom,hg19_pos,hg38_pos,top_trait,consort,ref,alt\n")
            hits_df = pd.read_csv(input.hit_list, dtype=str)

            for index, row in hits_df.iterrows():
                hg38_chrom = row['chrom']
                hg38_pos = row['pos']
                
                lifted, hg19_chrom, hg19_pos = lift_from_hg38_to_hg19(hg38_chrom, hg38_pos)
                assert lifted, f"[error] all the hit variants should be liftable = {hg38_chrom}:{hg38_pos}"

                top_trait = hg19_variant_to_trait[f"{hg19_chrom}:{hg19_pos}"]
                consort = "|".join(trait_to_consortium[top_trait])

                out_fd.write(f"{index},{hg19_chrom},{hg19_pos},{hg38_pos},{top_trait},{consort},{row['ref']},{row['alt']}\n")

        print(f"[log] finished building the variant masterlist\n")

rule build_variant_masterlist_for_all_nonhit_loci_analysis:
    input:
        nonhit_list="preprocessed_snp_lists/stingseq_nonhits_snp_list.csv",
        trait_list="analyses/stingseq_datasets/all_gwas_loci/metadata/trait_to_consortium.csv"
    output:
        masterlist="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.nonhit_loci.csv"
    run:
        import os
        import sys
        import pandas as pd 
        from liftover import get_lifter
        
        converter = get_lifter('hg19', 'hg38', one_based=True)
        ###########################################################
        # section 0: load in various dictionaries
        #            - trait-to-consortiums
        #            - hg19_variant-to-trait
        ###########################################################

        # load in a list of traits and what consortiums has those gwas
        trait_df = pd.read_csv(input.trait_list)
        trait_to_consortium = {}

        for _, row in trait_df.iterrows():
            trait_to_consortium[row['trait']] = row['consortium_list'].split("|")
        
        print(f"\n[log] loaded in list of traits (n = {len(trait_to_consortium)} unique traits)")

        # load in variant metadata to map variant to its gwas trait
        variants_df = pd.read_csv(os.path.join(data_dir, "stingseq_data/variant_list/john_table3f_v3.csv"), dtype=str)
        hg19_variant_to_trait = {}

        # cases where the variant file has a gwas trait that doesn't match table S1A
        special_cases = {"EO": "EOS",
                         "LYMPH_P": "LYM_P",
                         "NEUT": "NEU",
                         "EO_P": "EOS_P",
                         "LYMPH": "LYM",
                         "NEUT_P": "NEU_P",
                         "MONO": "MON",
                         "RBC_W": "RDW", # verified this one by hand
                         "BASO_P": "BAS_P"}

        # creates the variant to trait map
        for _, row in variants_df.iterrows():
            top_trait = row['Top GWAS fine-mapped']
            assert top_trait in special_cases or top_trait in trait_to_consortium, f"[error] unexpected trait found = {top_trait}"
            if top_trait in special_cases:
                top_trait = special_cases[top_trait]
            
            hg19_variant_to_trait[row['SNP Coordinates (hg19)']] = top_trait
            
        print("[log] finished building the dictionaries\n")

        ###########################################################
        # section 1: write variant masterlist
        ###########################################################

        def lift_from_hg38_to_hg19(hg38_chrom, hg38_pos):
            converter = get_lifter('hg38', 'hg19', one_based=True)
            pos_lifted = converter[f"chr{hg38_chrom}"][int(hg38_pos)]

            if len(pos_lifted) == 0:
                return False, "", ""
            else:
                hg19_chrom = pos_lifted[0][0].replace("chr", "")
                hg19_pos = int(pos_lifted[0][1])
                assert hg19_chrom == hg38_chrom

                return True, hg19_chrom, hg19_pos

        with open(output.masterlist, "w") as out_fd:
            out_fd.write("variant_id,chrom,hg19_pos,hg38_pos,top_trait,consort,ref,alt\n")
            nonhits_df = pd.read_csv(input.nonhit_list, dtype=str)

            for index, row in nonhits_df.iterrows():
                hg38_chrom = row['chrom']
                hg38_pos = row['pos']
                
                lifted, hg19_chrom, hg19_pos = lift_from_hg38_to_hg19(hg38_chrom, hg38_pos)
                assert lifted, f"[error] all the nonhit variants should be liftable = {hg38_chrom}:{hg38_pos}"

                top_trait = hg19_variant_to_trait[f"{hg19_chrom}:{hg19_pos}"]
                consort = "|".join(trait_to_consortium[top_trait])

                out_fd.write(f"{index},{hg19_chrom},{hg19_pos},{hg38_pos},{top_trait},{consort},{row['ref']},{row['alt']}\n")

        print(f"[log] finished building the variant masterlist (n = {nonhits_df.shape[0]})\n")


#######################################################
# section 1: extract variants from either hit or 
#            nonhit loci
#######################################################

rule extract_variants_from_ukbb_loci:
    input:
        masterlist="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.{loci_type}.csv"
    output:
        extracted_variants="analyses/stingseq_datasets/all_gwas_loci/sumstats/{loci_type}/variant{id}-{trait}-ukbb-variants.hg19.csv"
    shell:
        """
        # set the consortium variable
        consort="ukbb"

        # look through masterlist and verify these variables make sense
        awk -F',' -v id="{wildcards.id}" -v trait="{wildcards.trait}" -v consort="$consort" '
            NR==1 {{ next }}  # skip header
            $1 == id {{
                found_id = 1
                if ($5 != trait) {{
                    print "[error] trait mismatch for id {wildcards.id}"
                    exit 2
                }}
                n = split($6, a, "|")
                for (i = 1; i <= n; i++) {{
                    if (a[i] == consort) {{
                        exit 0
                    }}
                }}
                print "[error] consort '$consort' not found for id {wildcards.id}"
                exit 3
            }}
            END {{
                if (!found_id) {{
                    print "[error] id {wildcards.id} not found in masterlist"
                    exit 1
                }}
            }}
            ' {input.masterlist}
        
        # extract hit snp to use for querying
        hitsnp_chrom=$(awk -F',' '$1 == {wildcards.id}' {input.masterlist} | cut -d, -f2)
        hitsnp_pos=$(awk -F',' '$1 == {wildcards.id}' {input.masterlist} | cut -d, -f3)
        
        snp_id="${{hitsnp_chrom}}:${{hitsnp_pos}}"

        # build path to sumstats file to query
        gwas_sumstats_path="{data_dir}/stingseq_data/gwas/ukbb/{wildcards.trait}.tsv"

        if [ ! -f "$gwas_sumstats_path" ]; then
            echo "[error] GWAS sumstats file $gwas_sumstats_path does not exist"
            exit 1
        fi

        # extract variants from ukbb sumstats
        echo "chrom,pos,allele1,allele2,pval" > {output.extracted_variants}
        cat $gwas_sumstats_path | awk -F'\t' -v query=$snp_id -v window=500000 -v pthresh=5e-08 '
                                    BEGIN {{
                                        split(query, q, ":")
                                        qchr = q[1]
                                        qpos = q[2]
                                    }}
                                    NR == 1 {{ next }}   # skip header
                                    {{
                                        # Parse first column: chrom:pos:allele1:allele2
                                        split($1, a, ":")
                                        chr = a[1]
                                        pos = a[2]
                                        allele1 = a[3]
                                        allele2 = a[4]

                                        pval = $NF

                                        if (chr == qchr &&
                                            pos >= qpos - window &&
                                            pos <= qpos + window) {{
                                            print chr "," pos "," allele1 "," allele2 "," pval
                                        }}
                                    }}
                                    ' >> {output.extracted_variants}
        """

rule extract_variants_from_bcc_loci:
    input:
        masterlist="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.{loci_type}.csv"
    output:
        extracted_variants="analyses/stingseq_datasets/all_gwas_loci/sumstats/{loci_type}/variant{id}-{trait}-bcc-variants.hg19.csv"
    shell:
        """
        # set the consortium variable
        consort="bcc"

        # look through masterlist and verify these variables make sense
        awk -F',' -v id="{wildcards.id}" -v trait="{wildcards.trait}" -v consort="$consort" '
            NR==1 {{ next }}  # skip header
            $1 == id {{
                found_id = 1
                if ($5 != trait) {{
                    print "[error] trait mismatch for id {wildcards.id}"
                    exit 2
                }}
                n = split($6, a, "|")
                for (i = 1; i <= n; i++) {{
                    if (a[i] == consort) {{
                        exit 0
                    }}
                }}
                print "[error] consort '$consort' not found for id {wildcards.id}"
                exit 3
            }}
            END {{
                if (!found_id) {{
                    print "[error] id {wildcards.id} not found in masterlist"
                    exit 1
                }}
            }}
            ' {input.masterlist}
        
        # extract hit snp to use for querying
        hitsnp_chrom=$(awk -F',' '$1 == {wildcards.id}' {input.masterlist} | cut -d, -f2)
        hitsnp_pos=$(awk -F',' '$1 == {wildcards.id}' {input.masterlist} | cut -d, -f3)
        
        snp_id="${{hitsnp_chrom}}:${{hitsnp_pos}}"

        # build path to sumstats file to query
        gwas_sumstats_path="{data_dir}/stingseq_data/gwas/bcc/{wildcards.trait}.tsv.gz"

        if [ ! -f "$gwas_sumstats_path" ]; then
            echo "[error] GWAS sumstats file $gwas_sumstats_path does not exist"
            exit 1
        fi

        # extract variants from ukbb sumstats
        echo "chrom,pos,allele1,allele2,pval" > {output.extracted_variants}
        zcat $gwas_sumstats_path | awk -F'\t' -v query=$snp_id -v window=500000 -v pthresh=5e-08 '
                                    BEGIN {{
                                        split(query, q, ":")
                                        qchr = q[1]
                                        qpos = q[2]
                                    }}
                                    NR == 1 {{ next }}   # skip header
                                    {{
                                        # Parse first column: chrom:pos_allele1_allele2
                                        split($1, a, ":")
                                        chr = a[1]

                                        split(a[2], b, "_")
                                        pos = b[1]

                                        allele1 = $2
                                        allele2 = $3
                                        pval = $10

                                        if (chr == qchr &&
                                            pos >= qpos - window &&
                                            pos <= qpos + window) {{
                                            print chr "," pos "," allele1 "," allele2 "," pval
                                        }}
                                    }}
                                    ' >> {output.extracted_variants}
        """

#######################################################
# section 2: output generation rules to generate all 
#            needed hg19 sumstats for hit loci
#######################################################

def generate_list_of_all_possible_hg19_sumstats_for_hit_loci():
    """ returns a list of all possible sumstats we can extract """
    import os
    import pandas as pd

    masterlist_path = "analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.hit_loci.csv"
    if not os.path.exists(masterlist_path):
        return []
    
    input_files = []
    variant_masterlist_df = pd.read_csv(masterlist_path, dtype=str)
    
    for _, row in variant_masterlist_df.iterrows():
        for consort in row['consort'].split("|"):
            input_files.append(f"analyses/stingseq_datasets/all_gwas_loci/sumstats/hit_loci/variant{row['variant_id']}-{row['top_trait']}-{consort}-variants.hg19.csv")
    return input_files

rule generate_all_possible_stingseq_hit_hg19_sumstats:
    input:
        generate_list_of_all_possible_hg19_sumstats_for_hit_loci()

#######################################################
# section 3: output generation rules to generate all 
#            needed hg19 sumstats for nonhit loci
#######################################################

def generate_list_of_all_possible_hg19_sumstats_for_nonhit_loci():
    """ returns a list of all possible sumstats we can extract """
    import os
    import pandas as pd

    masterlist_path = "analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.nonhit_loci.csv"
    if not os.path.exists(masterlist_path):
        return []
    
    input_files = []
    variant_masterlist_df = pd.read_csv(masterlist_path, dtype=str)
    
    for _, row in variant_masterlist_df.iterrows():
        for consort in row['consort'].split("|"):
            input_files.append(f"analyses/stingseq_datasets/all_gwas_loci/sumstats/nonhit_loci/variant{row['variant_id']}-{row['top_trait']}-{consort}-variants.hg19.csv")
    return input_files

rule generate_all_possible_stingseq_nonhit_hg19_sumstats:
    input:
        generate_list_of_all_possible_hg19_sumstats_for_nonhit_loci()

#######################################################
# section 4: analyze the number of variants in each 
#            locus using different p-value cutoffs
#######################################################

rule analyze_variant_count_in_each_stingseq_sumstats:
    """
    this rule is meant to quickly analyze the sumstats
    after extracting them from the full gwas to see
    what would be a good cutoff to use for p-value

    that is why it called `*prior_to_finalizing.csv` 
    since I want to use it for exploratory analysis.
    """
    input:
        masterlist="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.{loci_type}.csv",
        hg19_hit_sumstats=generate_list_of_all_possible_hg19_sumstats_for_hit_loci(),
        hg19_nonhit_sumstats=generate_list_of_all_possible_hg19_sumstats_for_nonhit_loci()
    output:
        masterlist="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_counts_per_locus.{loci_type}.prior_to_finalizing.csv"
    run:
        import os 
        import sys
        import pandas as pd

        # open output file
        out_fd = open(output.masterlist, "w")
        out_fd.write("variant_id,chrom,pos,consort,pval_threshold,hit_snp_present,num_variants\n")
        
        # go through all variants and print out number of variants at each threshold
        variant_masterlist_df = pd.read_csv(input.masterlist, dtype=str)
        for index, row in variant_masterlist_df.iterrows():
            variant_id = row['variant_id']
            chrom = row['chrom']
            pos = row['hg19_pos']
            trait = row['top_trait']
            consort_list = row['consort'].split("|")

            for consort in consort_list:
                sumstats_path = f"analyses/stingseq_datasets/all_gwas_loci/sumstats/{wildcards.loci_type}/variant{variant_id}-{trait}-{consort}-variants.hg19.csv"
                sumstats_df = pd.read_csv(sumstats_path, dtype=str)

                for threshold in [5e-8, 1e-6]:
                    sumstats_df_subset = sumstats_df[sumstats_df['pval'].astype(float) <= threshold]
                    hit_snp_present = ((sumstats_df_subset['chrom'].astype(str) == str(chrom)) & (sumstats_df_subset['pos'].astype(str) == str(pos))).any()

                    out_fd.write(f"{variant_id},{chrom},{pos},{consort},{threshold},{hit_snp_present},{sumstats_df_subset.shape[0]}\n")

#######################################################
# section 5: postprocess the hg19 sumstats in order to
#            prepare them for deep learning scoring
#            - lift to hg38
#            - validate alleles are correct (e.g.
#              ref should be correct)
#######################################################

rule postprocess_variants_from_raw_sumstats_for_scoring:
    input:
        masterlist="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.{loci_type}.csv",
        extracted_variants="analyses/stingseq_datasets/all_gwas_loci/sumstats/{loci_type}/variant{id}-{trait}-{consort}-variants.hg19.csv"
    output:
        processed_variants="analyses/stingseq_datasets/all_gwas_loci/sumstats/{loci_type}/variant{id}-{trait}-{consort}-variants.hg38.csv"
    run:
        import os
        import sys
        import pysam
        import pandas as pd
        from liftover import get_lifter

        hg38_fasta = pysam.FastaFile(hg38_genome)
        converter = get_lifter('hg19', 'hg38', one_based=True)

        ###################################################
        # section 0: define helper methods for code
        ###################################################
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

        ###################################################
        # section 1: load hit snps and set of all hit snps
        ###################################################
        hg19_hitsnp_chrom = ""; hg19_hitsnp_pos = ""; hg38_hit_set = set()
        variant_masterlist_df = pd.read_csv(input.masterlist, dtype=str)

        for _, row in variant_masterlist_df.iterrows():
            if row['variant_id'] == wildcards.id:
                hg19_hitsnp_chrom = f"{row['chrom']}"
                hg19_hitsnp_pos = f"{row['hg19_pos']}"
            hg38_hit_set.add(f"{row['chrom']}:{row['hg38_pos']}")
        
        if hg19_hitsnp_chrom == "" and hg19_hitsnp_pos == "":
            print(f"[error] the hit snp (id = {wildcards.id}) does not occur in masterlist")
            exit(1)

        print(f"\n[log] loaded in the hit snp (id = {wildcards.id}) and hit snp set\n")

        ###################################################
        # section 2: load the extracted variants from the
        #            the gwas. do some checks:
        #            - hit snp should be present
        #            - hit snp should be significant
        #            - lift all variants
        #            - check the allele1/allele2 should
        #              contain the correct ref allele
        ###################################################
        input_variants_df = pd.read_csv(input.extracted_variants, dtype=str)
        pval_threshold = 5e-08

        # open output file and write header
        out_fd = open(output.processed_variants, "w")

        # check that hit snp is present
        hit_snp_present = (((input_variants_df['chrom'] == hg19_hitsnp_chrom) & (input_variants_df['pos'] == hg19_hitsnp_pos)).any())
        if not hit_snp_present:
            print(f"[error] hit snp (chrom={hg19_hitsnp_chrom}, pos={hg19_hitsnp_pos}) not found in the extracted variants file\n")

            print(f"[log] writing dummy output file ...\n")
            out_fd.write("DUMMY_FILE-NO_HIT_SNP_IN_THIS_SUMSTATS\n")
            return

        # check if hit snp is significant or not
        hit_snp_row = input_variants_df[(input_variants_df['chrom'] == hg19_hitsnp_chrom) & (input_variants_df['pos'] == hg19_hitsnp_pos)]
        hit_snp_pval = float(hit_snp_row.iloc[0]['pval'])

        if hit_snp_pval >= pval_threshold:
            print(f"[error] not a significant p-value for the hit snp (chrom={hg19_hitsnp_chrom}, pos={hg19_hitsnp_pos}) in the extracted variants file\n")

            print(f"[log] writing dummy output file ...\n")
            out_fd.write("DUMMY_FILE-NO_SIGNIFICANT_PVAL_FOR_HIT_SNP\n")
            return
        
        included_set = set() # variants included
        processed_snp_set = set() # variants already looked at
        duplicate_snp_list = [] # snps that occur multiple times
        lift_fail_set = set() # when lifting to hg38, there is an issue
        indel_fail_set = set() # variant is an indel, which is possible for all tools
        incorrect_ref_set = set() # variant's allele is not correct
        not_gws_set = set() # variant is not genome-wide significant

        # iterate through hg19 variants, lift them and check their alleles
        out_fd.write("chrom,pos,ref,alt,pval,hit\n")
        for _, row in input_variants_df.iterrows():
            hg19_chrom = row['chrom']
            hg19_pos = row['pos']
            ref_allele = row['allele1']
            alt_allele = row['allele2']
            pval = row['pval']

            snp = f"{hg19_chrom}:{hg19_pos}"
            
            # check if we have already processed this snp
            if snp in processed_snp_set:
                duplicate_snp_list.append(snp)
                continue
            processed_snp_set.add(snp)

            # check if the pvalue is significant
            if float(pval) >= pval_threshold:
                not_gws_set.add(snp)
                continue

            # check the liftover
            success_lift, hg38_chrom, hg38_pos = lift_from_hg19_to_hg38(hg19_chrom, hg19_pos)
            if not success_lift:
                lift_fail_set.add(snp)
                continue

            # check it is snp and not an indel
            if len(ref_allele) > 1 or len(alt_allele) > 1:
                indel_fail_set.add(snp)
                continue
            
            # check the alleles are correct (could be wrong order from stingseq file)
            correct_ref_base = hg38_fasta.fetch(f"chr{hg38_chrom}", hg38_pos-1, hg38_pos).upper()
            if ref_allele != correct_ref_base and alt_allele != correct_ref_base:
                incorrect_ref_set.add(snp)
                continue
            
            # flip them if necessary
            if correct_ref_base == alt_allele:
                alt_allele = ref_allele
                ref_allele = correct_ref_base
            
            # write out the variant if it passed the checks
            hitsnp_status = 1 if f"{hg38_chrom}:{hg38_pos}" in hg38_hit_set else 0
            out_fd.write(f"{hg38_chrom},{hg38_pos},{ref_allele},{alt_allele},{pval},{hitsnp_status}\n")

            included_set.add(snp)

        print(f"[log] number of variants removed due to duplicates: {len(duplicate_snp_list)}")
        print(f"[log] number of variants removed due to liftover: {len(lift_fail_set - included_set)}")
        print(f"[log] number of variants removed due to indel: {len(indel_fail_set - included_set)}")
        print(f"[log] number of variants removed due to incorrect ref base: {len(incorrect_ref_set - included_set)}")
        print(f"[log] number of variants removed due to not a significant p-value: {len(not_gws_set- included_set)}\n")

        print(f"[log] finished writing out the processed sumstats file (n = {len(included_set)})\n")
        
        out_fd.close()

#######################################################
# section 6: output generation rules for hg38 sumstats
#######################################################

def generate_list_of_all_possible_hg38_sumstats():
    """ returns a list of all possible sumstats we can extract """
    import os
    import pandas as pd

    input_files = []

    for loci_type in ['hit_loci', 'nonhit_loci']:
        masterlist_path = f"analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.{loci_type}.csv"
        if not os.path.exists(masterlist_path):
            return []
        variant_masterlist_df = pd.read_csv(masterlist_path, dtype=str)
        
        for _, row in variant_masterlist_df.iterrows():
            for consort in row['consort'].split("|"):
                input_files.append(f"analyses/stingseq_datasets/all_gwas_loci/sumstats/{loci_type}/variant{row['variant_id']}-{row['top_trait']}-{consort}-variants.hg38.csv")
    return input_files

rule generate_all_possible_stingseq_hg38_sumstats:
    input:
        generate_list_of_all_possible_hg38_sumstats()

#######################################################
# section 7: generate some metadata files that take
#            note on which variant extractions where
#            successful and which were not.
#            - key to success being that the hit 
#              variant occurs in that locus
#######################################################

rule validate_successful_loci_extractions:
    input:
        all_possible_sumstats=generate_list_of_all_possible_hg38_sumstats(),
        variant_masterlist="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_masterlist.{loci_type}.csv"
    output:
        validated_masterlist="analyses/stingseq_datasets/all_gwas_loci/metadata/variant_counts_per_locus.{loci_type}.finalized.csv"
    run:
        import os
        import sys
        import pandas as pd 

        # open output file
        out_fd = open(output.validated_masterlist, "w")
        out_fd.write("variant_id,chrom,hg19_pos,hg38_pos,ref,alt,top_trait,consort,num_variants,path\n")

        # load in variant masterlist, and check if variant extraction was successful
        variant_df = pd.read_csv(input.variant_masterlist, dtype=str)
        for index, row in variant_df.iterrows():
            variant_id = row['variant_id']
            chrom = row['chrom']
            hg19_pos = row['hg19_pos']
            hg38_pos = row['hg38_pos']
            ref = row['ref']
            alt = row['alt']
            trait = row['top_trait']
            consort_list = row['consort'].split("|")

            for consort in consort_list:
                sumstats_dir = f"analyses/stingseq_datasets/all_gwas_loci/sumstats/{wildcards.loci_type}"
                curr_sumstats = f"{sumstats_dir}/variant{variant_id}-{trait}-{consort}-variants.hg38.csv"

                # count lines to see if it successful; 1 means it is just the dummy line
                with open(curr_sumstats, "rb") as in_fd:
                    num_lines = sum(1 for _ in in_fd)
                
                if num_lines == 1:
                    continue
                else:
                    num_variants = num_lines - 1
                    out_fd.write(f"{variant_id},{chrom},{hg19_pos},{hg38_pos},{ref},{alt},{trait},{consort},{num_variants},{curr_sumstats}\n")
                
        out_fd.close()

#######################################################
# section A1: rules to run the hg38 sumstats through
#            alphagenome
#######################################################

def generate_list_of_alphagenome_scores_for_hg38_sumstats():
    import os
    import pandas as pd
    
    input_files = []

    specific_consort = "ukbb"
    min_variants = 50
    max_num_loci_to_score = 60

    for loci_type in ['hit_loci', 'nonhit_loci']:
        # validate masterlist exists
        masterlist_path = f"analyses/stingseq_datasets/all_gwas_loci/metadata/variant_counts_per_locus.{loci_type}.finalized.csv"
        if not os.path.exists(masterlist_path):
            return []
        variant_masterlist_df = pd.read_csv(masterlist_path, dtype=str)

        # subset to relevant sumstats
        variant_masterlist_df = variant_masterlist_df[variant_masterlist_df['consort'].str.contains(specific_consort)]
        variant_masterlist_df = variant_masterlist_df[variant_masterlist_df['num_variants'].astype(int) >= min_variants]
        
        # take the top n loci (to reduce scoring time)
        variant_masterlist_df = variant_masterlist_df.head(max_num_loci_to_score)

        for _, row in variant_masterlist_df.iterrows():
            input_files.append(f"analyses/stingseq_datasets/all_gwas_loci/variant_scores/alphagenome/{loci_type}/variant{row['variant_id']}-{row['top_trait']}-{row['consort']}-variants.alphagenome_vep_scores.csv")
    return input_files

rule generate_all_possible_stingseq_hit_alphagenome_scores:
    input:
        generate_list_of_alphagenome_scores_for_hg38_sumstats()

#######################################################
# section A2: rules to run the hg38 sumstats through
#            borzoi
#######################################################

def generate_list_of_borzoi_scores_for_hg38_sumstats():
    import os
    import pandas as pd
    
    input_files = []

    specific_consort = "ukbb"
    min_variants = 50
    max_num_loci_to_score = 60

    for loci_type in ['hit_loci', 'nonhit_loci']:
        # validate masterlist exists
        masterlist_path = f"analyses/stingseq_datasets/all_gwas_loci/metadata/variant_counts_per_locus.{loci_type}.finalized.csv"
        if not os.path.exists(masterlist_path):
            return []
        variant_masterlist_df = pd.read_csv(masterlist_path, dtype=str)

        # subset to relevant sumstats
        variant_masterlist_df = variant_masterlist_df[variant_masterlist_df['consort'].str.contains(specific_consort)]
        variant_masterlist_df = variant_masterlist_df[variant_masterlist_df['num_variants'].astype(int) >= min_variants]
        
        # take the top n loci (to reduce scoring time)
        variant_masterlist_df = variant_masterlist_df.head(max_num_loci_to_score)

        for _, row in variant_masterlist_df.iterrows():
            input_files.append(f"analyses/stingseq_datasets/all_gwas_loci/variant_scores/borzoi/{loci_type}/variant{row['variant_id']}-{row['top_trait']}-{row['consort']}-variants.borzoi_vep_scores.csv")
    return input_files

rule generate_all_possible_stingseq_hit_borzoi_scores:
    input:
        generate_list_of_borzoi_scores_for_hg38_sumstats()

#######################################################
# section A3: rules to run the hg38 sumstats through
#            borzoi prime
#######################################################

def generate_list_of_borzoi_prime_scores_for_hg38_sumstats():
    import os
    import pandas as pd
    
    input_files = []

    specific_consort = "ukbb"
    min_variants = 50
    max_num_loci_to_score = 60

    for loci_type in ['hit_loci', 'nonhit_loci']:
        # validate masterlist exists
        masterlist_path = f"analyses/stingseq_datasets/all_gwas_loci/metadata/variant_counts_per_locus.{loci_type}.finalized.csv"
        if not os.path.exists(masterlist_path):
            return []
        variant_masterlist_df = pd.read_csv(masterlist_path, dtype=str)

        # subset to relevant sumstats
        variant_masterlist_df = variant_masterlist_df[variant_masterlist_df['consort'].str.contains(specific_consort)]
        variant_masterlist_df = variant_masterlist_df[variant_masterlist_df['num_variants'].astype(int) >= min_variants]
        
        # take the top n loci (to reduce scoring time)
        variant_masterlist_df = variant_masterlist_df.head(max_num_loci_to_score)

        for _, row in variant_masterlist_df.iterrows():
            input_files.append(f"analyses/stingseq_datasets/all_gwas_loci/variant_scores/borzoi_prime/{loci_type}/variant{row['variant_id']}-{row['top_trait']}-{row['consort']}-variants.borzoi_prime_vep_scores.csv")
    return input_files

rule generate_all_possible_stingseq_hit_borzoi_prime_scores:
    input:
        generate_list_of_borzoi_prime_scores_for_hg38_sumstats()

#######################################################
# section A4: rules to run the hg38 sumstats through
#            enformer
#######################################################

def generate_list_of_enformer_scores_for_hg38_sumstats():
    import os
    import pandas as pd
    
    input_files = []

    specific_consort = "ukbb"
    min_variants = 50
    max_num_loci_to_score = 60

    for loci_type in ['hit_loci', 'nonhit_loci']:
        # validate masterlist exists
        masterlist_path = f"analyses/stingseq_datasets/all_gwas_loci/metadata/variant_counts_per_locus.{loci_type}.finalized.csv"
        if not os.path.exists(masterlist_path):
            return []
        variant_masterlist_df = pd.read_csv(masterlist_path, dtype=str)

        # subset to relevant sumstats
        variant_masterlist_df = variant_masterlist_df[variant_masterlist_df['consort'].str.contains(specific_consort)]
        variant_masterlist_df = variant_masterlist_df[variant_masterlist_df['num_variants'].astype(int) >= min_variants]
        
        # take the top n loci (to reduce scoring time)
        variant_masterlist_df = variant_masterlist_df.head(max_num_loci_to_score)

        for _, row in variant_masterlist_df.iterrows():
            input_files.append(f"analyses/stingseq_datasets/all_gwas_loci/variant_scores/enformer/{loci_type}/variant{row['variant_id']}-{row['top_trait']}-{row['consort']}-variants.enformer_vep_scores.csv")
    return input_files

rule generate_all_possible_stingseq_hit_enformer_scores:
    input:
        generate_list_of_enformer_scores_for_hg38_sumstats()