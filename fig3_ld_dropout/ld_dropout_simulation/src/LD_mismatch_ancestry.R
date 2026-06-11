###Creating LD mismatch files (autoimmune)
###Taking sumstats and LD matrix files and checking the overlapping variants and percent of variants excluded
###Neha Saravanan
###Update: 6/12/25

#Libraries
library(dplyr)
library(data.table)

##base directory
base_dir <- "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed"
study_folders <- list.dirs(base_dir, recursive = FALSE)
study_folders <- study_folders[grep("preprocessed$", study_folders)]

#list initialize
results_list <- list()

#loop through each study folder
for (study_folder in study_folders) {
  
  ##study name
  study <- basename(study_folder)
  #0.5 Mb files
  ancestry_folder <- list.dirs(study_folder, recursive = FALSE, full.names = TRUE)
  ancestry_folder <- ancestry_folder[grepl("_0.5Mb$", ancestry_folder)]
  #extract ancestry
  ancestry <- sub("(.*)_0.5Mb$", "\\1", basename(ancestry_folder))

  sumstats_dir <- file.path(base_dir, study, paste0(ancestry, "_0.5Mb"), "ss")
  ld_dir <- file.path(base_dir, study, paste0(ancestry, "_0.5Mb"), "ld", paste0(study, "_0.5Mb_"))
  
  #sumstats files
  sumstats_files <- list.files(sumstats_dir, pattern = "*.txt", full.names = TRUE)
  
  #empty data frame with initialized columns 
  results_df <- data.frame(LOCUS = character(),
                           sumstats_rows = numeric(),
                           ld_rows = numeric(),
                           sumstats_not_in_ld_rows = numeric(),
                           fraction_sumstats_not_in_ld = numeric(),
                           ancestry = character(),
                           stringsAsFactors = FALSE)
  
  for (file in sumstats_files) {
    sumstats <- read.table(file, header = TRUE)
    
    ##get locus name from the file name
    file_name <- basename(file)
    LOCUS <- sub(".*_([1-9X][0-9]?\\.[0-9]+)\\.txt", "\\1", file_name)
    
    #get LD file based on LOCUS value
    ld_file <- paste0(ld_dir, LOCUS, "_matched.tsv.bgz")
    print(ld_file)
    
    if (file.exists(ld_file)) {
      ld <- read.table(ld_file, header = TRUE)
      required_columns <- c("position", "allele1", "allele2", "beta")
      if (!all(required_columns %in% names(sumstats))) {
        message("Missing required columns in sumstats for LOCUS ", LOCUS)
        next
      }
      if (!all(required_columns %in% names(ld))) {
        message("Missing required columns in ld for LOCUS ", LOCUS)
        next
      }
      
      #convert all to characters
      sumstats$allele1 <- as.character(sumstats$allele1)
      sumstats$allele2 <- as.character(sumstats$allele2)
      sumstats$beta <- as.character(sumstats$beta)
      ld$allele1 <- as.character(ld$allele1)
      ld$allele2 <- as.character(ld$allele2)
      ld$beta <- as.character(ld$beta)
      
      #populate sumstats_not_in_ld
      sumstats_not_in_ld <- sumstats %>%
        anti_join(ld, by = c("position", "allele1", "allele2", "beta"))
      
      #get nrow for from sumstats and ld and sumstats_not_in_ld
      sumstats_rows <- nrow(sumstats)
      ld_rows <- nrow(ld)
      sumstats_not_in_ld_rows <- nrow(sumstats_not_in_ld)
      
      #calulcate fraction
      fraction_sumstats_not_in_ld <- sumstats_not_in_ld_rows / sumstats_rows
      
      #store all in dataframe :))
      results_df <- rbind(results_df, data.frame(LOCUS = LOCUS,
                                                 sumstats_rows = sumstats_rows,
                                                 ld_rows = ld_rows,
                                                 sumstats_not_in_ld_rows = sumstats_not_in_ld_rows,
                                                 fraction_sumstats_not_in_ld = fraction_sumstats_not_in_ld,
                                                 ancestry = ancestry))
    } else {
      message(paste("LD file for LOCUS", LOCUS, "not found in study", study, "and ancestry", ancestry, ". Skipping this LOCUS."))
    }
  }
  
  #for the specific study save in results_list
  results_list[[study]] <- results_df
}


##Merging results 

#initialize dataframe
merged_df <- data.frame()

#go through results_list and populate the new dataframe by study 
for (study in names(results_list)) {

  result_df <- results_list[[study]]
  if (nrow(result_df) > 0) {
     #keep track of origin study 
    result_df$study <- study
    #merge df
    merged_df <- rbind(merged_df, result_df)
  } else {
    message(paste("Skipping empty dataframe for study:", study))
  }
}


#get ancestry averages
ancestry_averages <- merged_df %>%
  filter(ancestry != "EUR") %>%
  group_by(ancestry) %>%
  summarise(avg_fraction_sumstats_not_in_ld = mean(fraction_sumstats_not_in_ld, na.rm = TRUE))
print(ancestry_averages)


#unique study and ancestry
study_ancestry_averages <- merged_df %>%
  group_by(study, ancestry) %>%
  summarise(avg_fraction_sumstats_not_in_ld = mean(fraction_sumstats_not_in_ld, na.rm = TRUE))
print(study_ancestry_averages)

write.table(merged_df, 
            file = "/gpfs/commons/groups/sanjana_lab/nsaravanan/omar/all_loci_multi.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
write.table(study_ancestry_averages, 
            file = "/gpfs/commons/groups/sanjana_lab/nsaravanan/omar/study.ancestry.averages.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

saveRDS(results_list, file = "/gpfs/commons/groups/sanjana_lab/nsaravanan/omar/results_list.rds")


############ FOR SIGNIFICANT VARIANTS 
base_dir <- "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/nathan_completed"
study_folders <- list.dirs(base_dir, recursive = FALSE)
study_folders <- study_folders[grep("preprocessed$", study_folders)]

#initialize significant results_list
results_list_significant <- list()

# Loop through each study folder
for (study_folder in study_folders) {
  study <- basename(study_folder)
  ancestry_folder <- list.dirs(study_folder, recursive = FALSE, full.names = TRUE)
  ancestry_folder <- ancestry_folder[grepl("_0.5Mb$", ancestry_folder)]
  ancestry <- sub("(.*)_0.5Mb$", "\\1", basename(ancestry_folder))
  sumstats_dir <- file.path(base_dir, study, paste0(ancestry, "_0.5Mb"), "ss")
  ld_dir <- file.path(base_dir, study, paste0(ancestry, "_0.5Mb"), "ld", paste0(study, "_0.5Mb_"))
  sumstats_files <- list.files(sumstats_dir, pattern = "*.txt", full.names = TRUE)
  results_df <- data.frame(LOCUS = character(),
                           sumstats_rows = numeric(),
                           ld_rows = numeric(),
                           sumstats_not_in_ld_rows = numeric(),
                           fraction_sumstats_not_in_ld = numeric(),
                           ancestry = character(),
                           stringsAsFactors = FALSE)
  
  for (file in sumstats_files) {
    sumstats <- read.table(file, header = TRUE)
    ###SIGNIFICANCE FILTERING 
    sumstats <- sumstats[sumstats$pval <= 5e-8, ]
    file_name <- basename(file)
    LOCUS <- sub(".*_([1-9X][0-9]?\\.[0-9]+)\\.txt", "\\1", file_name)
    ld_file <- paste0(ld_dir, LOCUS, "_matched.tsv.bgz")
    print(ld_file)
    
    #ld 
    if (file.exists(ld_file)) {
      ld <- read.table(ld_file, header = TRUE)
      ###SIGNIFICANCE FILTERING
      ld <- ld %>% filter(ld$pval < 5e-8)
      required_columns <- c("position", "allele1", "allele2", "beta")
      if (!all(required_columns %in% names(sumstats))) {
        message("Missing required columns in sumstats for LOCUS ", LOCUS)
        next
      }
      if (!all(required_columns %in% names(ld))) {
        message("Missing required columns in ld for LOCUS ", LOCUS)
        next
      }
      
      #character conversion 
      sumstats$allele1 <- as.character(sumstats$allele1)
      sumstats$allele2 <- as.character(sumstats$allele2)
      sumstats$beta <- as.character(sumstats$beta)
      
      ld$allele1 <- as.character(ld$allele1)
      ld$allele2 <- as.character(ld$allele2)
      ld$beta <- as.character(ld$beta)
      sumstats_not_in_ld <- sumstats %>%
        anti_join(ld, by = c("position", "allele1", "allele2", "beta"))
      sumstats_rows <- nrow(sumstats)
      ld_rows <- nrow(ld)
      sumstats_not_in_ld_rows <- nrow(sumstats_not_in_ld)
      fraction_sumstats_not_in_ld <- sumstats_not_in_ld_rows / sumstats_rows
      results_df <- rbind(results_df, data.frame(LOCUS = LOCUS,
                                                 sumstats_rows = sumstats_rows,
                                                 ld_rows = ld_rows,
                                                 sumstats_not_in_ld_rows = sumstats_not_in_ld_rows,
                                                 fraction_sumstats_not_in_ld = fraction_sumstats_not_in_ld,
                                                 ancestry = ancestry))
    } else {
      message(paste("LD file for LOCUS", LOCUS, "not found in study", study, "and ancestry", ancestry, ". Skipping this LOCUS."))
    }
  }
  results_list_significant[[study]] <- results_df
}

##merge results
merged_df_1 <- data.frame()
for (study in names(results_list_significant)) {
  result_df <- results_list_significant[[study]]
  if (nrow(result_df) > 0) {
    result_df$study <- study
    merged_df_1 <- rbind(merged_df_1, result_df)
  } else {
    message(paste("Skipping empty dataframe for study:", study))
  }
}

ancestry_averages <- merged_df_1 %>%
  filter(ancestry != "EUR") %>%
  group_by(ancestry) %>%
  summarise(avg_fraction_sumstats_not_in_ld = mean(fraction_sumstats_not_in_ld, na.rm = TRUE))
print(ancestry_averages)

study_ancestry_averages_1 <- merged_df_1 %>%
  group_by(study, ancestry) %>%
  summarise(avg_fraction_sumstats_not_in_ld = mean(fraction_sumstats_not_in_ld, na.rm = TRUE))

write.table(merged_df_1, 
            file = "/gpfs/commons/groups/sanjana_lab/nsaravanan/omar/all_loci_multi_1.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
write.table(study_ancestry_averages_1, 
            file = "/gpfs/commons/groups/sanjana_lab/nsaravanan/omar/study.ancestry.averages_1.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

saveRDS(results_list_significant, file = "/gpfs/commons/groups/sanjana_lab/nsaravanan/omar/results_list_1.rds")

##Creating data frames summarizing loci for each study

merged_totals <- merged_df %>%
  group_by(study, ancestry) %>%
  summarise(
    sumstats_rows = sum(sumstats_rows, na.rm = TRUE),
    ld_rows = sum(ld_rows, na.rm = TRUE),
    sumstats_not_in_ld_rows = sum(sumstats_not_in_ld_rows, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    fraction_not_included = sumstats_not_in_ld_rows / sumstats_rows,
    percent_not_included = fraction_not_included * 100
  ) %>%
  arrange(desc(percent_not_included)) %>%
  slice(-c(1:5))

merged_totals_1 <- merged_df_1 %>%
  group_by(study, ancestry) %>%
  summarise(
    sumstats_rows = sum(sumstats_rows, na.rm = TRUE),
    ld_rows = sum(ld_rows, na.rm = TRUE),
    sumstats_not_in_ld_rows = sum(sumstats_not_in_ld_rows, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    fraction_not_included = sumstats_not_in_ld_rows / sumstats_rows,
    percent_not_included = fraction_not_included * 100
  ) %>%
  arrange(desc(percent_not_included)) %>%
  slice(-c(1:5))

write.table(merged_totals, file = "/gpfs/commons/groups/sanjana_lab/nsaravanan/omar/merged_totals.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(merged_totals_1, file = "/gpfs/commons/groups/sanjana_lab/nsaravanan/omar/merged_totals_1.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

##AGE check
merged_totals_age <- merged_totals

merged_totals_age$OG_IC <- ifelse(grepl("^2008|^2009|^2010|^2011|^2012", merged_totals_age$study), 1, 0)
merged_totals_age$part_IC <- ifelse(grepl("^2013|^2014|^2015|^2016", merged_totals_age$study), 1, 0)
merged_totals_age$arrays <- ifelse(grepl("^2017|^2018", merged_totals_age$study), 1, 0)
merged_totals_age$high_density_arrays <- ifelse(grepl("^2019|^2020|^2021|^2022|^2023|^2024|^2025", merged_totals_age$study), 1, 0)

write.table(merged_totals_age, file = "/gpfs/commons/groups/sanjana_lab/nsaravanan/omar/merged_totals_age.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
