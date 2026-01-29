######################################################
# Name: preprocess_stingseq_variants.smk
# Description: preprocesses variant list from the 
#              sting-seq paper
# Date: january 4th, 2025
######################################################

rule parse_stingseq_variants_to_lookup_in_dbsnp:
    output:
        variant_list="workspace/sting_seq/bcftools_input.tsv"
    run:
        import os
        import re
        import pandas as pd 

        ##############################################
        # section 0: load variant table from paper
        ##############################################
        variants_df = pd.read_csv(os.path.join(data_dir, "stingseq_data/variant_list/john_table3f_v3.csv"), 
                                  dtype=str)
        out_fd = open(output.variant_list, "w")

        ##################################################
        # section 1: extract variants that are hits and 
        #            identify variants only with rsids
        ##################################################

        # subset variants dataframe to ones with significant p-values
        hits_df = variants_df[variants_df["Q-value (500 kb)"].notna() & 
                             (variants_df["Q-value (500 kb)"] != "") &
                             (variants_df["Q-value (500 kb)"].astype(float) < 0.05)
                             ]
        
        # create a map from unique coordinates to SNP ids
        pos_to_ids = {} 
        raw_id_pattern = re.compile(r'^\d+:\d+_[ACGT]+_[ACGT]+$')
        rs_id_pattern = re.compile(f'rs\d+')

        for index, row in hits_df.iterrows():
            var_id = row['SNP Coordinates (hg19)']
            if var_id not in pos_to_ids:
                pos_to_ids[var_id] = {"rsid": "",
                                      "raw_id": ""}
            
            snp = str(row['SNP'])
            if raw_id_pattern.match(snp):
                pos_to_ids[var_id]['raw_id'] = snp
            elif rs_id_pattern.match(snp):
                pos_to_ids[var_id]['rsid'] = snp
            else:
                print(f"[error] unexpected snp structure = {snp}"); exit(1)
        
        # write out the hit snps with only rsid information
        for pos, data in pos_to_ids.items():
            if len(data['raw_id']) == 0:
                chrom = pos.split(":")[0]
                pos = pos.split(":")[1]
                out_fd.write(f"{chrom}\t{pos}\n")
        
        # convert the coordinate column to set 
        hit_coordinates_set = set(hits_df["SNP Coordinates (hg19)"])

        ##################################################
        # section 2: extract variants that are non-hits, 
        #            identify variants only with rsids and
        #            make sure it is not in the hit list
        ##################################################

        # subset variants dataframe to ones with significant p-values
        nonhits_df = variants_df[variants_df["Q-value (500 kb)"].notna() & 
                                (variants_df["Q-value (500 kb)"] != "") &
                                (variants_df["Q-value (500 kb)"].astype(float) >= 0.05)
                                ]
        
        # create a map from unique coordinates to SNP ids
        pos_to_ids = {} 
        for index, row in nonhits_df.iterrows():
            var_id = row['SNP Coordinates (hg19)']
            if var_id not in pos_to_ids:
                pos_to_ids[var_id] = {"rsid": "",
                                      "raw_id": ""}
            
            snp = str(row['SNP'])
            if raw_id_pattern.match(snp):
                pos_to_ids[var_id]['raw_id'] = snp
            elif rs_id_pattern.match(snp):
                pos_to_ids[var_id]['rsid'] = snp
            else:
                print(f"[error] unexpected snp structure = {snp}"); exit(1)
        
        # write out the nonhit snps with only rsid information and not a hit
        for pos, data in pos_to_ids.items():
            if len(data['raw_id']) == 0 and pos not in hit_coordinates_set:
                chrom = pos.split(":")[0]
                pos = pos.split(":")[1]
                out_fd.write(f"{chrom}\t{pos}\n")
        
        out_fd.close()

rule run_bcftools_to_extract_data_for_stingseq_snps:
    input:
        variant_list="workspace/sting_seq/bcftools_input.tsv"
    output:
        output_1="workspace/sting_seq/bcftools_output_1.vcf.gz",
        output_2="workspace/sting_seq/bcftools_output_2.vcf.gz",
        output_final="workspace/sting_seq/bcftools_output.vcf.gz",
    shell:
        """
        bcftools view -R {input.variant_list} {dbSNP_all} -Oz -o {output.output_1}
        bcftools view -R {input.variant_list} {dbSNP_common} -Oz -o {output.output_2}

        bcftools index {output.output_1}
        bcftools index {output.output_2}

        bcftools concat -a {output.output_1} {output.output_2} -Oz -o {output.output_final}
        bcftools index {output.output_final}
        """

rule generate_variant_lists_for_stingseq_snps:
    input:
        orig_variant_list="workspace/sting_seq/bcftools_input.tsv",
        dbsnp_results="workspace/sting_seq/bcftools_output.vcf.gz"
    output:
        stingseq_hits="preprocessed_snp_lists/stingseq_hits_snp_list.csv",
        stingseq_nonhits="preprocessed_snp_lists/stingseq_nonhits_snp_list.csv"
    run:
        import pysam
        import os
        import glob
        import pandas as pd
        from liftover import get_lifter
        
        converter = get_lifter('hg19', 'hg38', one_based=True)
        hg38_fasta = pysam.FastaFile(hg38_genome)

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


        ##################################################
        # section 1: load in bcftools output 
        ##################################################
        dbsnp_vcf = pysam.VariantFile(input.dbsnp_results)
        snp_to_alleles = {}

        valid_chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y']
        valid_alleles = ['A', 'C', 'G', 'T']

        # grab all the snps found in dbsnp and store in dictionary
        for record in dbsnp_vcf.fetch():
            curr_chr = str(record.chrom)
            assert curr_chr in valid_chromosomes, f"error: obtained invalid chromosome={curr_chr}"
            
            curr_pos = int(record.pos)
            curr_ref = str(record.ref)
            assert all(ch in valid_alleles for ch in curr_ref), f"error: obtained invalid ref allele={curr_ref}"

            curr_alleles = tuple(record.alts)
            assert all(ch in valid_alleles for x in curr_alleles for ch in x), f"error: obtained invalid alt alleles={curr_alleles}"

            var_id = f"{curr_chr}:{curr_pos}"
            if var_id not in snp_to_alleles:
                snp_to_alleles[var_id] = []
            
            for curr_alt in curr_alleles:
                if (curr_ref, curr_alt) not in snp_to_alleles[var_id]:
                    snp_to_alleles[var_id].append((curr_ref, curr_alt))
        
        print(f"\n[log] finished loading in the allele info from dbsnp\n")

        ##################################################
        # section 2: load in variants again from and merge
        #            in dbsnp information when needed
        ##################################################
        variants_df = pd.read_csv(os.path.join(data_dir, "stingseq_data/variant_list/john_table3f_v3.csv"), 
                                  dtype=str)

        # subset variants dataframe to ones with significant p-values
        hits_df = variants_df[variants_df["Q-value (500 kb)"].notna() & 
                             (variants_df["Q-value (500 kb)"] != "") &
                             (variants_df["Q-value (500 kb)"].astype(float) >= 0.00) &
                             (variants_df["Q-value (500 kb)"].astype(float) < 0.05)
                             ]
        
        # print starting number and keep track of losses
        print(f"[log] aftering filtering on q-value, we have {hits_df['SNP Coordinates (hg19)'].nunique()} hits")
        
        included_set = set() # variants included
        lift_fail_set = set() # when lifting to hg38, there is an issue
        indel_fail_set = set() # variant is an indel, which is possible for all tools
        incorrect_ref_set = set() # variant's allele is not correct

        # create a map from unique coordinates to SNP ids
        pos_to_ids = {} 
        raw_id_pattern = re.compile(r'^\d+:\d+_[ACGT]+_[ACGT]+$')
        rs_id_pattern = re.compile(f'rs\d+')

        for index, row in hits_df.iterrows():
            var_id = row['SNP Coordinates (hg19)']
            if var_id not in pos_to_ids:
                pos_to_ids[var_id] = {"rsid": "",
                                      "raw_id": ""}
            
            snp = str(row['SNP'])
            if raw_id_pattern.match(snp):
                pos_to_ids[var_id]['raw_id'] = snp
            elif rs_id_pattern.match(snp):
                pos_to_ids[var_id]['rsid'] = snp
            else:
                print(f"[error] unexpected snp structure = {snp}"); exit(1)
        
        ###############################################################
        # hard-coded alleles: for some of the sting-seq hits, there 
        #                     multiple alleles and the only way to 
        #                     figure out the correct one is to track
        #                     it down from the sumstats. I manually 
        #                     curated these to make sure we have the
        #                     correct allele
        ###############################################################
        special_hits_cases = {"11:67804807": ('G', 'A'),
                              "12:6441623": ('A', 'G'),
                              "12:6503786": ('T', 'C'),
                              "13:114147674": ('A', 'G'),
                              "18:9191860": ('A', 'C'),
                              "19:45413233": ('G', 'T'),
                              "5:134722833": ('C', 'A'),
                              "6:109625797": ('A', 'G'),
                              "8:11720227": ('G', 'C')}

        # write out the hits list (rely on raw_id, then use dbsnp info)
        with open(output.stingseq_hits, "w") as out_fd:
            out_fd.write("chrom,pos,ref,alt\n")

            for snp, data in pos_to_ids.items():
                if len(data['rsid']) and len(data['raw_id']) == 0:
                    assert snp in snp_to_alleles, f"[error] no dbsnp information for this snp = {pos}"
                    hg19_chrom = snp.split(":")[0]
                    hg19_pos = snp.split(":")[1]

                    # check the liftover
                    success_lift, hg38_chrom, hg38_pos = lift_from_hg19_to_hg38(hg19_chrom, hg19_pos)
                    if not success_lift:
                        lift_fail_set.add(snp)
                        continue
                    
                    # handle special cases
                    allele_list = snp_to_alleles[snp].copy()
                    if snp in special_hits_cases:
                        allele_list = [special_hits_cases[snp]]

                    # go through all alleles (filter to only snps)
                    for ref_allele, alt_allele in allele_list:
                        # check if it is a snp and not an indel
                        if len(ref_allele) > 1 or len(alt_allele) > 1:
                            indel_fail_set.add(snp)
                            continue

                        # check the alleles are correct (could be wrong order from stingseq file)
                        correct_ref_base = hg38_fasta.fetch(f"chr{hg38_chrom}", hg38_pos-1, hg38_pos).upper()
                        if ref_allele != correct_ref_base and alt_allele != correct_ref_base:
                            incorrect_ref_set.add(snp)
                            continue
                            
                        if correct_ref_base == alt_allele:
                            alt_allele = ref_allele
                            ref_allele = correct_ref_base

                        included_set.add(snp)
                        out_fd.write(f"{hg38_chrom},{hg38_pos},{ref_allele},{alt_allele}\n")
                        
                        # TO-DO: here I am breaking after printing 1 allele, but I should choose a specific allele
                        break

                elif len(data['raw_id']):
                    hg19_chrom = data['raw_id'].split(":")[0]
                    hg19_pos = data['raw_id'].split(":")[1].split("_")[0]

                    alt_allele = data['raw_id'].split("_")[1]
                    ref_allele = data['raw_id'].split("_")[2]

                    assert hg19_chrom == snp.split(":")[0] and hg19_pos == snp.split(":")[1], f"[error] incongruence between {snp} and {data['raw_id']}"

                    # check the liftover
                    success_lift, hg38_chrom, hg38_pos = lift_from_hg19_to_hg38(hg19_chrom, hg19_pos)
                    if not success_lift:
                        lift_fail_set.add(snp)
                        continue
                    
                    # check the alleles are correct (could be wrong order from stingseq file)
                    correct_ref_base = hg38_fasta.fetch(f"chr{hg38_chrom}", hg38_pos-1, hg38_pos).upper()
                    if ref_allele != correct_ref_base and alt_allele != correct_ref_base:
                        incorrect_ref_set.add(snp)
                        continue
                    
                    if correct_ref_base == alt_allele:
                        alt_allele = ref_allele
                        ref_allele = correct_ref_base

                    # check it is snp and not an indel
                    if len(ref_allele) > 1 or len(alt_allele) > 1:
                        indel_fail_set.add(snp)
                        continue
                    
                    included_set.add(snp)
                    out_fd.write(f"{hg38_chrom},{hg38_pos},{ref_allele},{alt_allele}\n")
        
        print(f"\n[log] number of variants removed due to liftover: {len(lift_fail_set - included_set)}")
        print(f"[log] number of variants removed due to indel: {len(indel_fail_set - included_set)}")
        print(f"[log] number of variants removed due to incorrect ref base: {len(incorrect_ref_set - included_set)}\n")
        
        print(f"[log] finished writing out the variant list for stingseq hits\n")

        # convert the coordinate column to set 
        hit_coordinates_set = set(hits_df["SNP Coordinates (hg19)"])

        ##################################################
        # section 3: merge in dbsnp info for nonhits and
        #            and then write them out
        ##################################################

        # subset variants dataframe to ones with significant p-values
        nonhits_df = variants_df[variants_df["Q-value (500 kb)"].notna() & 
                                (variants_df["Q-value (500 kb)"] != "") &
                                (variants_df["Q-value (500 kb)"].astype(float) >= 0.05)
                                ]
        print(f"[log] aftering filtering on q-value, we have {nonhits_df['SNP Coordinates (hg19)'].nunique()} non-hits")

        included_set = set() # variants included
        overlap_hit_set = set() # variants already in the hit set
        lift_fail_set = set() # when lifting to hg38, there is an issue
        indel_fail_set = set() # variant is an indel, which is possible for all tools
        incorrect_ref_set = set() # variant's allele is not correct

        # create a map from unique coordinates to SNP ids
        pos_to_ids = {} 
        for index, row in nonhits_df.iterrows():
            var_id = row['SNP Coordinates (hg19)']
            if var_id not in pos_to_ids:
                pos_to_ids[var_id] = {"rsid": "",
                                      "raw_id": ""}
            
            snp = str(row['SNP'])
            if raw_id_pattern.match(snp):
                pos_to_ids[var_id]['raw_id'] = snp
            elif rs_id_pattern.match(snp):
                pos_to_ids[var_id]['rsid'] = snp
            else:
                print(f"[error] unexpected snp structure = {snp}"); exit(1)

        # write out the nonhits list (rely on raw_id, then use dbsnp info)
        with open(output.stingseq_nonhits, "w") as out_fd:
            out_fd.write("chrom,pos,ref,alt\n")

            for snp, data in pos_to_ids.items():
                # ignore snps that do have significant associations
                if snp in hit_coordinates_set:
                    overlap_hit_set.add(snp)
                    continue

                # check if we should use rsid or raw_id
                if len(data['rsid']) and len(data['raw_id']) == 0:
                    assert snp in snp_to_alleles, f"[error] no dbsnp information for this snp = {pos}"
                    hg19_chrom = snp.split(":")[0]
                    hg19_pos = snp.split(":")[1]

                    # check the liftover
                    success_lift, hg38_chrom, hg38_pos = lift_from_hg19_to_hg38(hg19_chrom, hg19_pos)
                    if not success_lift:
                        lift_fail_set.add(snp)
                        continue

                    # go through all alleles (filter to only snps)
                    for ref_allele, alt_allele in snp_to_alleles[snp]:
                        # check if it is a snp and not an indel
                        if len(ref_allele) > 1 or len(alt_allele) > 1:
                            indel_fail_set.add(snp)
                            continue

                        # check the alleles are correct (could be wrong order from stingseq file)
                        correct_ref_base = hg38_fasta.fetch(f"chr{hg38_chrom}", hg38_pos-1, hg38_pos).upper()
                        if ref_allele != correct_ref_base and alt_allele != correct_ref_base:
                            incorrect_ref_set.add(snp)
                            continue
                            
                        if correct_ref_base == alt_allele:
                            alt_allele = ref_allele
                            ref_allele = correct_ref_base

                        included_set.add(snp)
                        out_fd.write(f"{hg38_chrom},{hg38_pos},{ref_allele},{alt_allele}\n")

                        # TO-DO: here I am breaking after printing 1 allele, but I should choose a specific allele
                        break

                elif len(data['raw_id']):
                    hg19_chrom = data['raw_id'].split(":")[0]
                    hg19_pos = data['raw_id'].split(":")[1].split("_")[0]

                    alt_allele = data['raw_id'].split("_")[1]
                    ref_allele = data['raw_id'].split("_")[2]

                    assert hg19_chrom == snp.split(":")[0] and hg19_pos == snp.split(":")[1], f"[error] incongruence between {snp} and {data['raw_id']}"

                    # check the liftover
                    success_lift, hg38_chrom, hg38_pos = lift_from_hg19_to_hg38(hg19_chrom, hg19_pos)
                    if not success_lift:
                        lift_fail_set.add(snp)
                        continue

                    # check the alleles are correct (could be wrong order from stingseq file)
                    correct_ref_base = hg38_fasta.fetch(f"chr{hg38_chrom}", hg38_pos-1, hg38_pos).upper()
                    if ref_allele != correct_ref_base and alt_allele != correct_ref_base:
                        incorrect_ref_set.add(snp)
                        continue
                    
                    if correct_ref_base == alt_allele:
                        alt_allele = ref_allele
                        ref_allele = correct_ref_base
                    
                    # check it is snp and not an indel
                    if len(ref_allele) > 1 or len(alt_allele) > 1:
                        indel_fail_set.add(snp)
                        continue

                    included_set.add(snp)
                    out_fd.write(f"{hg38_chrom},{hg38_pos},{ref_allele},{alt_allele}\n")

        print(f"\n[log] number of variants removed due overlapping hit set: {len(overlap_hit_set - included_set)}")
        print(f"[log] number of variants removed due to liftover: {len(lift_fail_set - included_set)}")
        print(f"[log] number of variants removed due to indel: {len(indel_fail_set - included_set)}")
        print(f"[log] number of variants removed due to incorrect ref base: {len(incorrect_ref_set - included_set)}\n")

        print(f"[log] finished writing out the variant list for stingseq non-hits\n")