# library(rmarkdown)
library(dplyr)
library(stringr)
# library(knitr)

# library(digest)

# Genetic data ----
data_path <- "../../scicore_mirror/data/acmguru_data_20250107"
dir.exists(data_path)
list.files(path = data_path, full.names = TRUE)

# Small tables ----
file_list_small <- list.files(path = data_path, pattern = "ACMGuru_singlecase_df_report_main_text_v1_momic_samplecount_180_chr_", full.names = TRUE)

# Large tables ----
file_list <- list.files(path = data_path, pattern = "ACMGuru_singlecase_df_report_v1_momic_samplecount_180_chr_", full.names = TRUE)

# Read all RDS files and bind them into a single dataframe
dt_selected_sub <- file_list %>%
# dt_selected_sub <- file_list_small %>%
  lapply(readRDS) %>%
  bind_rows()

# Splitting the string based on "NGS" and creating a new column
dt_selected_sub <- dt_selected_sub %>%
  mutate(sample.id_meta = str_extract(sample.id, "NGS.*"),
         sample.id = str_remove(sample.id, "NGS.*"))

# Get the best variant ID
dt_selected_sub <- dt_selected_sub %>% mutate(HGVSpc = coalesce(HGVSp, HGVSc))

dt_selected_sub <- dt_selected_sub %>% filter(AF.x < 0.3)


# tests
df_mthfr <-  dt_selected_sub|> filter(SYMBOL == "MTHFR")
df_mmachc <-  dt_selected_sub|> filter(SYMBOL == "MMACHC")
df_tcn2 <-  dt_selected_sub|> filter(SYMBOL == "TCN2")
df_mmachc$HGVSpc |> unique()
df_mthfr$HGVSpc |> unique()
df_tcn2$HGVSpc |> unique()

# # expected vs small
# MMACHC	2 variants (55 cases)	 - 4 found
# MTHFR (39 cases)	38 variants	- 22 found
# TCN2 (1 case) - 6 found
# 
# # expected versus autodetected
# MMACHC	2 variants (55 cases)	 - 5 found variants
# MTHFR (39 cases)	38 variants	- 23 found variants
# TCN2 (1 case) - 6 found variants

# Knwon
# ENST00000215838.8:c.172del;ENST00000405742.7:c.172del;ENST00000407817.3:c.172del;ENST00000423350.1:n.501del;ENST00000450638.5:c.97del;ENST00000471659.2:n.343del;ENST00000698263.1:c.172del;ENST00000698264.1:n.343del;ENST00000698265.1:c.172del;ENST00000698266.1:c.172del;ENST00000698267.1:c.172del;ENST00000698268.1:c.172del;ENST00000698269.1:c.172del;ENST00000698270.1:c.172del;ENST00000698271.1:c.172del;ENST00000698272.1:c.172del;ENST00000698273.1:c.163del	ENSP00000215838.3:p.Leu58TyrfsTer28;ENSP00000385914.3:p.Leu58TyrfsTer28;ENSP00000384914.3:p.Leu58TyrfsTer28;ENSP00000394184.2:p.Leu33TyrfsTer28;ENSP00000513635.1:p.Leu58TyrfsTer28;ENSP00000513636.1:p.Leu58TyrfsTer28;ENSP00000513637.1:p.Leu58TyrfsTer28;ENSP00000513638.1:p.Leu58TyrfsTer28;ENSP00000513639.1:p.Leu58TyrfsTer28;ENSP00000513640.1:p.Leu58TyrfsTer28;ENSP00000513641.1:p.Leu58TyrfsTer28;ENSP00000513642.1:p.Leu58TyrfsTer28;ENSP00000513643.1:p.Leu58TyrfsTer28;ENSP00000513644.1:p.Leu55TyrfsTer28 = missing note non-canonical

# ENST00000401061.9:c.271dup = Yes
# ENST00000401061.9:c.394C>T = Yes

# filters
# dt_selected_sub <- dt_selected_sub |> filter(gnomAD_AF < 0.4)
# dt_selected_sub <- dt_selected_sub |> filter(MAX_AF < 0.4)

# Splitting the string based on "NGS" and creating a new column
dt_selected_sub <- dt_selected_sub %>%
  mutate(sample.id_meta = str_extract(sample.id, "NGS.*"),
         sample.id = str_remove(sample.id, "NGS.*"))

# # To print reports in order of highest rank guru score
# # Step 1: Group by 'sample.id'
# # Step 2: Calculate the max 'ACMG_total_score' for each 'sample.id'
# # Step 3: Arrange the dataframe based on the max score
# # Calculate max score for each sample and rank them
# dt_selected_sub_ranked_samples <- dt_selected_sub %>%
#   group_by(sample.id) %>%
#   summarise(max_score_in_case = max(ACMG_total_score, na.rm = TRUE)) %>%
#   arrange(desc(max_score_in_case)) %>%
#   mutate(case_rank_in_cohort = row_number()) %>%  # Add a rank column
#   ungroup()
# 
# # Join back to get the other columns
# dt_selected_sub_ranked <- dt_selected_sub_ranked_samples %>%
#   left_join(dt_selected_sub, by = "sample.id")
# 
# dt_selected_sub <- dt_selected_sub_ranked
# rm(dt_selected_sub_ranked, dt_selected_sub_ranked_samples)
# 


# To print reports in order of highest rank guru score
# Step 1: Group by 'sample.id'
# Step 2: Calculate the max 'ACMG_total_score' for each 'sample.id'
# Step 3: Arrange the dataframe based on the max score
# Calculate max score for each sample and rank them

print("This step causes some known hits to drop in MTHFR, MMACHC since they don't have the highest score.")

# dt_selected_sub_ranked_samples <- dt_selected_sub %>%
#   group_by(sample.id) %>%
#   summarise(max_score_in_case = max(ACMG_total_score, na.rm = TRUE)) %>%
#   arrange(desc(max_score_in_case)) %>%
#   mutate(case_rank_in_cohort = row_number()) %>%  # Add a rank column
#   ungroup()
# 
# # Join back to get the other columns
# dt_selected_sub_ranked <- dt_selected_sub_ranked_samples %>%
#   left_join(dt_selected_sub, by = "sample.id")



print("Rule: Special genes that we know cause THIS disease.")

# # Group by sample.id and consider special handling for MMACHC or MTHFR
# 
# dt_selected_sub_ranked <- dt_selected_sub %>%
#   # Include a column to identify special genes
#   mutate(is_special_gene = SYMBOL %in% c("MMACHC", "MTHFR")) %>%
#   # Group by sample.id
#   group_by(sample.id) %>%
#   # Calculate max score within each group
#   mutate(max_score_in_case = max(ACMG_total_score, na.rm = TRUE)) %>%
#   ungroup() %>%
#   # Filter to keep either the max score or special gene rows
#   filter(ACMG_total_score == max_score_in_case | is_special_gene) %>%
#   # Arrange by descending max score
#   arrange(desc(max_score_in_case)) %>%
#   # Create a rank column based on the arranged data
#   mutate(case_rank_in_cohort = row_number())

library(dplyr)

library(dplyr)

library(dplyr)

library(dplyr)

library(dplyr)

# First, rank variants within each case
dt_selected_sub <- dt_selected_sub %>%
  group_by(sample.id) %>%
  # Calculate the maximum score within each group to identify the top priority variant
  mutate(max_score_in_case = max(ACMG_total_score, na.rm = TRUE)) %>%
  ungroup()

  
library(dplyr)

dt_selected_sub_ranked <- dt_selected_sub %>%
  mutate(is_special_gene = SYMBOL %in% c("MMACHC", "MTHFR")) %>%
  group_by(sample.id) %>%
  mutate(variant_rank_in_case = dense_rank(desc(is_special_gene) + desc(ACMG_total_score))) %>%
  ungroup()   %>%
dplyr::select("sample.id", 
              # "case_rank_in_cohort", 
              "variant_rank_in_case",
              "ACMG_total_score", "ACMG_count", "source", "SYMBOL", "HGVSp", "HGVSc", 
              "Consequence", "IMPACT", "genotype", "gnomAD_AF", "comp_het_flag", everything())

# dt_selected_sub_ranked2 <- dt_selected_sub_ranked %>%
#   mutate(is_special_gene = SYMBOL %in% c("MMACHC", "MTHFR")) %>%
#   group_by(sample.id) %>%
#   mutate(case_rank_in_cohort = dense_rank(desc(is_special_gene) + desc(ACMG_total_score))) %>%
#   ungroup()   %>%
#   dplyr::select("sample.id", 
#                 "case_rank_in_cohort",
#                 "variant_rank_in_case",
#                 "ACMG_total_score", "ACMG_count", "source", "SYMBOL", "HGVSp", "HGVSc", 
#                 "Consequence", "IMPACT", "genotype", "gnomAD_AF", "comp_het_flag", everything())
# 
















# 
# 
# 
# library(dplyr)
# 
# # Filter to top variants where variant_rank_in_case == 1
# top_variants <- dt_selected_sub_ranked %>%
#   filter(variant_rank_in_case == 1)
# 
# library(dplyr)
# 
# # Assuming 'top_variants' contains the correct 'case_rank_in_cohort' calculations
# top_variants <- top_variants %>%
#   mutate(is_special_gene = SYMBOL %in% c("MMACHC", "MTHFR")) %>%
#   mutate(case_rank_in_cohort = dense_rank(desc(is_special_gene) + desc(ACMG_total_score))) %>%
#   dplyr::select(sample.id, case_rank_in_cohort_top = case_rank_in_cohort)  # Rename for clarity
# 
# # Join the calculated case_rank_in_cohort back to the original dataset
# dt_selected_sub_ranked2 <- dt_selected_sub_ranked %>%
#   left_join(top_variants, by = "sample.id") %>%
#   mutate(case_rank_in_cohort = coalesce(case_rank_in_cohort_top, case_rank_in_cohort)) %>%
#   dplyr::select(-case_rank_in_cohort_top)  # Remove the temporary rank column
# 
# 
# 
# dt_selected_sub_ranked2 <- merge(dt_selected_sub_ranked, top_variants, by = "sample.id", all.x = T)
#   
#   
# # View the updated data frame to confirm changes
# print(dt_selected_sub_ranked2)
# 
# 
# 
# 
# # %>%
#   dplyr::select("sample.id", "case_rank_in_cohort", "variant_rank_in_case",
#          "ACMG_total_score", "ACMG_count", "source", "SYMBOL", "HGVSp", "HGVSc", 
#          "Consequence", "IMPACT", "genotype", "gnomAD_AF", "comp_het_flag", everything())
# 
# 
# 
# # Optionally, filter to keep only the top ranked variants per sample.id if needed later
# top_ranked_variants <- dt_selected_sub_ranked %>%
#   filter(variant_rank_in_case == 1)
# 
# 
# d_test <- dt_selected_sub_ranked |> filter(sample.id == "HCY040")
# d_test <- dt_selected_sub_ranked |> filter(sample.id == "HCY140")
# 
# 
# # Optionally, filter to keep only the top ranked variants per sample.id if needed later
# top_ranked_variants <- dt_selected_sub_ranked %>%
#   filter(case_rank_in_cohort == 1)



# Now, dt_selected_sub_ranked_samples includes rows where either they are a max scorer or are a special gene.

# Clean up: Optionally remove the temporary columns if they are no longer needed
# dt_selected_sub_ranked_samples <- dt_selected_sub_ranked_samples %>%
  # select(-c(max_score_in_case, is_special_gene))

# Now, dt_selected_sub_ranked_samples will keep variants of MMACHC or MTHFR even if they aren't the highest scoring.





# # Calculate within-case rank for each variant
# dt_selected_sub_ranked <- dt_selected_sub_ranked %>%
#   group_by(sample.id) %>%
#   mutate(variant_rank_in_case = rank(-ACMG_total_score, ties.method = "first")) %>%
#   ungroup()
# 
# # Calculate within-case rank for each variant where equal scores are ranked the same
# dt_selected_sub_ranked <- dt_selected_sub_ranked %>%
#   group_by(sample.id) %>%
#   mutate(variant_rank_in_case = min_rank(-ACMG_total_score)) %>%
#   ungroup()
# 
# 
# dt_selected_sub <- dt_selected_sub_ranked
# rm(dt_selected_sub_ranked, dt_selected_sub_ranked_samples)
# 
# dt_selected_sub <- dt_selected_sub |>
#   dplyr::select("sample.id", 
#          "max_score_in_case",
#          "case_rank_in_cohort",
#          "variant_rank_in_case",
#        "ACMG_total_score",
#        # "ACMG_score", 
#        "ACMG_count",
#        "source",
#        # "ACMG_highest", 
#        "SYMBOL",
#        # "chr", #"rownames",
#        "HGVSp", 
#        "HGVSc",
#        "Consequence", 
#        # "cohort_pheno",
#        # "frequency_in_1", 
#        "IMPACT",
#        "genotype", 
#        # "Inheritance",
#        # "CLIN_SIG", 
#        "gnomAD_AF",
#        "comp_het_flag",   
# 
#        everything()
#        # "SymbolCount"
# ) 

# dt_selected_sub_filter06 <- dt_selected_sub |> filter(ACMG_total_score > 5) 
# dt_selected_sub_filter10 <- dt_selected_sub |> filter(ACMG_total_score > 9) 
dt_selected_sub_rank1 <- dt_selected_sub_ranked |> filter(variant_rank_in_case < 2) 
dt_selected_sub_rank2 <- dt_selected_sub_ranked |> filter(variant_rank_in_case < 3) 

data.table::fwrite(dt_selected_sub_ranked, file="../../scicore_mirror/data/acmguru_data_20250107/ACMGuru_singlecase_df_report_v1_momic_samplecount_180_all_chr_filterNone.csv", sep = ",")

# data.table::fwrite(dt_selected_sub_filter06, file="../../scicore_mirror/data/acmguru_data_20250107/ACMGuru_singlecase_df_report_v1_momic_samplecount_180_all_chr_filter06.csv", sep = ",")

# data.table::fwrite(dt_selected_sub_filter10, file="../../scicore_mirror/data/acmguru_data_20250107/ACMGuru_singlecase_df_report_v1_momic_samplecount_180_all_chr_filter10.csv", sep = ",")

data.table::fwrite(dt_selected_sub_rank1, file="../../scicore_mirror/data/acmguru_data_20250107/ACMGuru_singlecase_df_report_v1_momic_samplecount_180_all_chr_rank1.csv", sep = ",")

data.table::fwrite(dt_selected_sub_rank2, file="../../scicore_mirror/data/acmguru_data_20250107/ACMGuru_singlecase_df_report_v1_momic_samplecount_180_all_chr_rank2.csv", sep = ",")

