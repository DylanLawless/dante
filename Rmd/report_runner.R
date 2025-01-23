library(rmarkdown)
library(dplyr)
library(dplyr)
library(knitr)
library(stringr)
library(digest)
# kableExtra.latex.load_packages=FALSE
# library(kableExtra)

case_count_test <- 1
print(paste("Number of cases to run:", case_count_test))

# Import main report table.
# Get top score per patient, group by gene, 
# Get discussion
print("Look at MaxScore and descide how many candidates to report.")
# prepare data ----
# single file test
dt_selected_sub <-  readRDS("../../scicore_mirror/data/acmguru_data_20250107/ACMGuru_singlecase_df_report_main_text_v1_momic_samplecount_180_chr_1.Rds")

# Guru figure data ----
pdf_path <- "../../scicore_mirror/data/acmguru_data_20250107/pdf"
dir.exists(pdf_path)
list.files(path = pdf_path, full.names = TRUE)
pdf_list <-list.files(path = pdf_path, full.names = TRUE)

# Genetic data ----
data_path <- "../../scicore_mirror/data/acmguru_data_20250107"
dir.exists(data_path)
list.files(path = data_path, full.names = TRUE)
file_list <- list.files(path = data_path, pattern = "ACMGuru_singlecase_df_report_main_text_v1_momic_samplecount_180_chr_", full.names = TRUE)

# Read all RDS files and bind them into a single dataframe
dt_selected_sub <- file_list %>%
  lapply(readRDS) %>%
  bind_rows()

# To print reports in order of highest rank guru score
# Step 1: Group by 'sample.id'
# Step 2: Calculate the max 'ACMG_total_score' for each 'sample.id'
# Step 3: Arrange the dataframe based on the max score
# Calculate max score for each sample and rank them
dt_selected_sub_ranked_samples <- dt_selected_sub %>%
  group_by(sample.id) %>%
  summarise(max_score = max(ACMG_total_score, na.rm = TRUE)) %>%
  arrange(desc(max_score)) %>%
  mutate(rank = row_number()) %>%  # Add a rank column
  ungroup()

# Join back to get the other columns
dt_selected_sub_ranked <- dt_selected_sub_ranked_samples %>%
  left_join(dt_selected_sub, by = "sample.id")

dt_selected_sub <- dt_selected_sub_ranked
rm(dt_selected_sub_ranked, dt_selected_sub_ranked_samples)

names(dt_selected_sub)

# ACMGuru_singlecase_df_report_v1_momic_samplecount_180_chr_1
# ACMGuru_singlecase_df_report_main_text_v1_momic_samplecount_180_chr_1

# Splitting the string based on "NGS" and creating a new column
dt_selected_sub <- dt_selected_sub %>%
  mutate(sample.id_meta = str_extract(sample.id, "NGS.*"),
         sample.id = str_remove(sample.id, "NGS.*"))

dt_selected_sub$sample.id %>% unique()
dt_selected_sub$sample.id_meta %>% unique()

maindataset <- dt_selected_sub %>%
  mutate(
    HGVSc = str_replace(HGVSc, ":", ": "),
    HGVSp = str_replace(HGVSp, ":", ": "),
    # Consequence = str_replace_all(Consequence, "_", "- ")
    Consequence = str_replace_all(Consequence, "_", " "),
    Consequence = str_replace_all(Consequence, "&", " and ") # required to escape the & in latex tables
  ) %>%
  # mutate(across(c(HGVSc, HGVSp), ~str_wrap(.x, 10))) %>%
  # mutate(across(everything(), ~gsub("\n", "<br> ", .x))) %>% # Replace \n with LaTeX \\ for line breaks
  group_by(sample.id) %>%
  # Compute the max ACMG_total_score for each group
  mutate(MaxScore = max(ACMG_total_score)) %>%
  # Filter to keep only the rows where the score equals the maximum score
  filter(ACMG_total_score == MaxScore) %>% 
  # Ungroup to perform operations not confined by the group
  ungroup() %>%
  # Count the number of unique SYMBOLs per sample.id
  add_count(sample.id, SYMBOL, name = "SymbolCount") %>%
  # Filter to keep only those groups where there's more than one SYMBOL
  group_by(sample.id) %>%
  # filter(any(SymbolCount > 1)) %>% 
  ungroup() %>%
  select("sample.id", 
         "ACMG_total_score",
         # "ACMG_score", 
         "ACMG_count",
         # "ACMG_highest", 
         "SYMBOL",
         # "chr", #"rownames",
         "HGVSp", 
         "HGVSc",
         "Consequence", 
         # "cohort_pheno",
         # "frequency_in_1", 
         "IMPACT",
         "genotype", 
         # "Inheritance",
         # "CLIN_SIG", 
         "gnomAD_AF",
         "comp_het_flag",   
         "MaxScore", 
         "rank"
         # "SymbolCount"
  ) 

# replace chars for latex ----
replace_chars_for_latex <- function(x) {
  stringr::str_replace_all(x, c("&" = "and", "_" = " ", "\\+" = "$+$"))
}

maindataset <- maindataset %>%
  mutate_if(is.character, replace_chars_for_latex)

maindataset <- maindataset %>%
  arrange(rank)

# colnames(df)[colnames(df) == 'oldName'] <- 'newName'
# colnames(maindataset)[colnames(maindataset) == 'ACMG_total_score'] <- 'ACMG total score'
# colnames(maindataset)[colnames(maindataset) == 'ACMG_count'] <- 'ACMG count'
# colnames(maindataset)[colnames(maindataset) == 'comp_het_flag'] <- 'comp het flag'

# View the results
head(maindataset)

library(dplyr)
library(knitr)
# library(kableExtra)

# tabular report variables ----

# run report ----
print("Note that a new 'df' would be sourced in report.Rmd so do not use this variable without checking")


# Deduplicate to get unique sample.ids for report generation
unique_ids <- unique(maindataset$sample.id)

# TEST limit to 2 samples ----
unique_ids <- head(unique_ids, case_count_test)  # limited test set

# Loop through each unique sample.id
for(sample_id in unique_ids) {
  
  print(paste("Now running report: ", sample_id))
  print(paste("..............................."))
  
  # Filter data for the current sample_id
  current_data <- maindataset[maindataset$sample.id == sample_id, ]
  
  # Serialize the dataframe to get a consistent format for checksum calculation
  serialized_data <- serialize(current_data, NULL)
  data_checksum <- digest(serialized_data, algo = "md5")
  
  print(paste("Data checksum for sample ID", sample_id, ":", data_checksum))
  
  # Assuming maindataset is already loaded and available
  current_data <- current_data %>%
    filter(sample.id == sample.id) %>%
    slice(1:5) %>%
    # mutate(across(everything(), ~str_wrap(.x, 10))) %>%
    # t() %>%
    as.data.frame()
  
  # Extract the rank of the current sample
  current_rank <- current_data$rank |> unique()
  print(paste("rank:", current_rank))
  # Define the output file path with rank in the filename
  report_path <- paste0("../reports/report_priority_", current_rank, "_sample_", sample_id , ".pdf")
  # report_path <- paste0("../reports/report_", sample_id, ".pdf")
  
  render("report.Rmd", output_file = report_path,
         params = list(
           singlecase_ID = sample_id,
           data_checksum = data_checksum,
           pdf_list = pdf_list
         ),
         quiet = TRUE)
  
  log_file <-  paste0("./report_priority_", current_rank, "_sample_", sample_id , ".log")
  file.remove(log_file)  # Removes the log file
  
  print(paste("Done."))
}
