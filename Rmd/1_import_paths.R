# General paths ----
pca_path <- "../data/PCA/" 
case_out <- "../output/" 

# Guru figure data ----
pdf_path <- "../data/guru/pdf"
dir.exists(pdf_path)
list.files(path = pdf_path, full.names = TRUE)
pdf_list <-list.files(path = pdf_path, full.names = TRUE)

# Genetic data ----
data_path <- "../data/guru"
dir.exists(data_path)
list.files(path = data_path, full.names = TRUE)
file_list <- list.files(path = data_path, pattern = "ACMGuru_singlecase_df_raw_v1_momic_samplecount_180_chr_", full.names = TRUE)

# cleaning processing ----
path_clean_unprocessed <- "../data/clean_unprocessed/"
path_clean_unprocessed_impact <- "../data/clean_unprocessed_impact/"
path_preprocessed <- "../data/preprocessed/"
