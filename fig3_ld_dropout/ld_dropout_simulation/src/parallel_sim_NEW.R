library(data.table)
library(susieR)
library(Matrix)
library(dplyr)

.libPaths("/gpfs/commons/home/nsaravanan/R/x86_64-pc-linux-gnu-library/4.4")
library(GWASBrewer)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript sim_v3.R <locus_idx> <h2_val> <N_val> <pi_val>")
}

#-----------------------------------------------
# Parse command-line arguments
#-----------------------------------------------
locus_idx <- as.integer(args[1])
h2_val <- as.numeric(args[2])
N_val <- as.numeric(args[3])
pi_val <- as.numeric(args[4])

print(paste("Locus Index:", locus_idx))
print(paste("Heritability:", h2_val))
print(paste("Sample Size:", N_val))
print(paste("Pi:", pi_val))

#-----------------------------------------------
# Set paths and functions
#-----------------------------------------------
set.seed(42)

generate_af <- function(n) rbeta(n, 1, 5)

ld_matrix_paths <- c(
  "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed/2020_32296059_AST_EUR_preprocessed/EUR_0.5Mb/ld/UKBB_LDmat_2020_32296059_AST_EUR_preprocessed_0.5Mb_2.72054837.bgz",
  "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed/2020_32296059_AST_EUR_preprocessed/EUR_0.5Mb/ld/UKBB_LDmat_2020_32296059_AST_EUR_preprocessed_0.5Mb_10.6093139.bgz",
  "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed/2020_32296059_AST_EUR_preprocessed/EUR_0.5Mb/ld/UKBB_LDmat_2020_32296059_AST_EUR_preprocessed_0.5Mb_1.8505058.bgz",
  "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed/2024_PanUKBB_UC_EUR_preprocessed/EUR_0.5Mb/ld/UKBB_LDmat_2024_PanUKBB_UC_EUR_preprocessed_0.5Mb_21.40464924.bgz",
  "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed/2024_PanUKBB_PSO_EUR_preprocessed/EUR_0.5Mb/ld/UKBB_LDmat_2024_PanUKBB_PSO_EUR_preprocessed_0.5Mb_6.26184733.bgz",
  "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed/2021_33536424_SLE_EUR-EAS_preprocessed/EUR-EAS_0.5Mb/ld/UKBB_LDmat_2021_33536424_SLE_EUR-EAS_preprocessed_0.5Mb_5.150457485.bgz",
  "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed/2022_35896530_SD_EUR_preprocessed/EUR_0.5Mb/ld/UKBB_LDmat_2022_35896530_SD_EUR_preprocessed_0.5Mb_8.10828909.bgz",
  "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed/2024_PanUKBB_CD_EUR_preprocessed/EUR_0.5Mb/ld/UKBB_LDmat_2024_PanUKBB_CD_EUR_preprocessed_0.5Mb_16.50660964.bgz",
  "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed/2024_PanUKBB_SLE_EUR_preprocessed/EUR_0.5Mb/ld/UKBB_LDmat_2024_PanUKBB_SLE_EUR_preprocessed_0.5Mb_21.37058961.bgz",
  "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed/2024_PanUKBB_PSO_EUR_preprocessed/EUR_0.5Mb/ld/UKBB_LDmat_2024_PanUKBB_PSO_EUR_preprocessed_0.5Mb_11.128401358.bgz"
)

#-----------------------------------------------
# Load and prepare LD matrix
#-----------------------------------------------
ld_matrix_path <- ld_matrix_paths[locus_idx]
ld_dt <- fread(ld_matrix_path)
ld_matrix <- as.matrix(ld_dt)
ld_matrix[lower.tri(ld_matrix)] <- t(ld_matrix)[lower.tri(ld_matrix)]
nearPD_mat <- Matrix::nearPD(ld_matrix, conv.tol=1e-10)
R <- as.matrix(nearPD_mat$mat)
af_vector <- generate_af(nrow(R))
num_variants <- nrow(R)

#-----------------------------------------------
# Run 10 GWAS simulations for this h2, N, pi
#-----------------------------------------------
all_results <- list()
dropout_percents <- c(0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80)

for (sim_rep in 1:10) {
  set.seed(1000 + sim_rep)  # ensure reproducibility
  
  sim_dat <- sim_mv(G = 1,
                    J = num_variants,
                    N = N_val,
                    h2 = h2_val,
                    pi = pi_val/(num_variants),
                    af = af_vector,
                    pi_exact = TRUE,
                    est_s = TRUE,
                    R_LD = list(R))
  
  beta_hat <- sim_dat$beta_hat[, 1]
  se_hat <- sim_dat$se_beta_hat[, 1]
  z_scores <- beta_hat / se_hat
  causal_idx <- which(sim_dat$direct_SNP_effects_joint != 0)
  
  for (dropout in dropout_percents) {
    num_drop <- ceiling((dropout / 100) * num_variants)
    
    for (rep_dropout in 1:3) {
      drop_idx <- sort(sample(seq_len(num_variants), size=num_drop, replace=FALSE))
      keep_idx <- setdiff(seq_len(num_variants), drop_idx)
      
      susie_fit <- susie_rss(
        bhat = beta_hat[keep_idx],
        shat = se_hat[keep_idx],
        R = R[keep_idx, keep_idx],
        n = N_val,
        L = 10
      )
      
      PIPs <- susie_fit$pip
      full_result <- data.frame(
        variant_index = keep_idx,
        replicate = sim_rep,
        dropout_percent = dropout,
        dropout_replicate = rep_dropout,
        causal = ifelse(keep_idx %in% causal_idx, "Y", "N"),
        PIP = PIPs
      )
      
      all_results[[length(all_results) + 1]] <- full_result
    }
  }
}

#-----------------------------------------------
# Write merged CSV with all runs
#-----------------------------------------------
output_df <- do.call(rbind, all_results)

outdir <- "/gpfs/commons/groups/sanjana_lab/nsaravanan/sim_v6"
dir.create(outdir, showWarnings = FALSE)

outfile <- file.path(outdir, paste0("PIPs_locus", locus_idx, "_h2_", h2_val, "_N_", N_val, "_pi_", round(pi_val, 6), ".csv"))
fwrite(output_df, outfile)
print(paste("[log] Finished. Output saved to:", outfile))