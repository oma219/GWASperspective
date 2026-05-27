# Identify colocalization examples for paper:

In this README, I describe the steps I took to identify interesting loci to use for the paper. In case, I want to redo it or change it slighly I know what files I looked at exactly.

## Step 1: Identify loci from the overall file

I have this file `data/signal_list/unique_signals.step_10.csv` that contains details for all unique signals in our dataset. The relevant features/columns in this file are `singleton`, `colocalized`, `quant_only`, and `num_gwas`. 

Using these features, I can get a list of signals where it colocalized across multiple traits and then it is not only quantiative traits but also we want as many GWAS as possible. Here is a list of the columns in this file for reference:

```txt
1 chromosome
2 start
3 end
4 region_size
5 singleton
6 colocalized
7 quant_only
8 num_gwas
9 gwas_with_leadsnps
10 num_gwas_with_leadsnps
11 gwas_without_leadsnps
12 num_gwas_without_leadsnps
13 num_disease_traits
14 num_protein_traits
15 num_blood
16 variants
17 num_variants
18 num_fm_ld_variants
19 num_coloc_variants
20 num_hyprcoloc_variants
21 num_recovery_variants
```

I ran this following command where I wanted something with these conditions:

- Not a singleton signal
- Colocalized signal across different GWAS
- Not only contains quantitative traits
- No nominated variants from finemapping
- At least 1 variant nominated from using hyprcoloc

This command should show two columns where the first one corresponds to number of gwas colocalized and the second column is the number of disease traits. So essentially I was looking at those rows and trying to identify signals with a large number of GWAS and large number of unique diseases. 

```sh
awk -F, '$5 == "False" && $6 == "True" && $7 == "False" && $18 == 0 && $20 > 0' unique_signals.step_10.csv | sort -t, -k8,8nr | cut -d, -f8,13 | less -N
```

### Example 1: Signal with 16 GWAS and 9 diseases

Command: `awk -F, '$5 == "False" && $6 == "True" && $7 == "False" && $18 == 0 && $20 > 0' unique_signals.step_10.csv | sort -t, -k8,8nr | head -n8 | tail -n1 | awk -F, '{for(i=1;i<=NF; i++){print i, $i;}}'`

Fields:

```txt
1 4
2 39791916
3 40807564
4 1015648
5 False
6 True
7 False
8 16
9 4.40307564|2020_32514122_GD_EAS/4.40307564|2021_34594039EAS_GD_EAS/4.40319646|2024_38982041_AITD_EUR/4.40307564|2021_34594039_GD_EUR-EAS/4.40291916|2020_32581359_AITD_EUR/4.40307564|2016_27723757_VIT_EUR
10 6
11 2019_31604244_MS_EUR|2020_32005708_T1D_EUR|2022_36333501EAS_RA_EAS|2020_33310728_RA_EUR-EAS|2021_34594039_RA_EUR-EAS|2022_35470158_RA_EUR|2021_34012112_T1D_EUR|2020_33106285_JIA_EUR|2023_37156999_CD_EUR-EAS|2016_27723758_SIGD_EUR
12 10
13 9
14 0
15 0
16 4.40307564
17 1
18 0
19 1
20 1
21 0
```

This example is interesting for various reasons. We see a high number of diseases colocalizing and we see a good number of GWAS with and without leadsnps showing the power of colocalization. In addition, we see that no SNPs from this locus was nominated from finemapping but we identify a SNP with hyprcoloc.

For referring to the hyprcoloc outputs, I looked up this region using the start and end coordinates to identify which region id this was that I tried to colocalize by looking in this file: `/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/colocalization/v2_scripts/hyprcoloc_region_boundaries_with_signals.txt`

And I identified that this region is `563` and here is the full line of information from that file:

```txt
563	4	40191916	40407564	1049:1050:1051:1052:1053:1054
```


### Example 2: Signal with 10 GWAS and 8 diseases

Command: `awk -F, '$5 == "False" && $6 == "True" && $7 == "False" && $18 == 0 && $20 > 0' unique_signals.step_10.csv | sort -t, -k8,8nr | head -n60 | tail -n1 | awk -F, '{for(i=1;i<=NF; i++){print i, $i;}}'`

Fields:
```txt
1 7
2 4936753
3 5936753
4 1000000
5 False
6 True
7 False
8 10
9 7.5436753|2021_34594039_IBD_EUR-EAS/7.5436753|2023_37156999_IBD_EUR-EAS/7.5436753|2023_37156999_CD_EUR-EAS/7.5436753|2023_37156999_UC_EUR-EAS
10 4
11 2021_34594039_PV_EUR-EAS|2021_34594039_T1D_EUR-EAS|2021_34594039_HT_EUR-EAS|2020_32581359_AITD_EUR|2024_38982041_AITD_EUR|2021_34594039_GD_EUR-EAS
12 6
13 8
14 0
15 0
16 7.5436753|7.5162926
17 2
18 0
19 1
20 2
21 0
```

Again, good to see a large number of GWAS (10) and a large number of distinct diseases (8). And we have a decent number of diseases that do not have lead snps so again it shows the power of colocalization. And again we see now SNPs nominated by finemapping across any of these GWAS but with colocalization we identify two variants. 

For referring to the hyprcoloc outputs, I looked up this region using the start and end coordinates to identify which region id this was that I tried to colocalize by looking in this file: `/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/colocalization/v2_scripts/hyprcoloc_region_boundaries_with_signals.txt`

And I identified that this region is `677` and here is the full line of information from that file:

```txt
677	7	5336753	5536753	2493:2494:2495:2496
```

## Step 2: Copy over sumstats files for each example

### Example 1: Signal with 16 GWAS and 9 diseases

```sh
cd data/example_sumstats/example_1

cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2020_32514122_GD_EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2021_34594039EAS_GD_EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2024_38982041_AITD_EUR_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2021_34594039_GD_EUR-EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2020_32581359_AITD_EUR_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2016_27723757_VIT_EUR_preprocessed.sorted.tsv .

cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2019_31604244_MS_EUR_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2020_32005708_T1D_EUR_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2022_36333501EAS_RA_EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2020_33310728_RA_EUR-EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2021_34594039_RA_EUR-EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2022_35470158_RA_EUR_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2021_34012112_T1D_EUR_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2020_33106285_JIA_EUR_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2023_37156999_CD_EUR-EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2016_27723758_SIGD_EUR_preprocessed.sorted.tsv .
```

### Example 2: Signal with 10 GWAS and 8 diseases

```sh
cd data/example_sumstats/example_2

cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2021_34594039_IBD_EUR-EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2023_37156999_IBD_EUR-EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2023_37156999_CD_EUR-EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2023_37156999_UC_EUR-EAS_preprocessed.sorted.tsv .

cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2021_34594039_PV_EUR-EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2021_34594039_T1D_EUR-EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2021_34594039_HT_EUR-EAS_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2020_32581359_AITD_EUR_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2024_38982041_AITD_EUR_preprocessed.sorted.tsv .
cp /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/preprocessed_v2/2021_34594039_GD_EUR-EAS_preprocessed.sorted.tsv .
```

## Step 3: Extract relevant sumstats for plotting

### Example 1: Signal with 16 GWAS and 9 diseases

```sh
# set these variables in the bash script
INPUT_DIR="/gpfs/commons/groups/sanjana_lab/oahmed/perspective_coloc_example/data/example_sumstats/example_1"

REGION_CHR="4"
REGION_START=39791916
REGION_END=40807564

OUTPUT_PATH="/gpfs/commons/groups/sanjana_lab/oahmed/perspective_coloc_example/results/example_1/region_sumstats.csv"
> $OUTPUT_PATH
```

```sh
./extract_relevant_sumstats.sh
```

### Example 2: Signal with 10 GWAS and 8 diseases

```sh
# set these variables in the bash script
INPUT_DIR="/gpfs/commons/groups/sanjana_lab/oahmed/perspective_coloc_example/data/example_sumstats/example_2"

REGION_CHR="7"
REGION_START=4936753
REGION_END=5936753

OUTPUT_PATH="/gpfs/commons/groups/sanjana_lab/oahmed/perspective_coloc_example/results/example_2/example2_region_sumstats.csv"
> $OUTPUT_PATH
```