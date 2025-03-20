library(rmarkdown)
library(dplyr)
library(knitr)
library(stringr)
library(digest)
library(ggplot2)
library(ggrepel)

# minimal set ---- 
col_names <- c("sample.id","ACMG_total_score",
               "SYMBOL","HGVSp",
               "HGVSc","max_score",
               "rank",
               "rownames","ACMG_count",
               "Consequence","IMPACT",
               "genotype", "gnomAD_AF",
               "comp_het_flag","genotype_name",
               "ACMG_criteria","MaxScore"
               )

# custom -----
prior_score_weight <- 20 # known disease panel
case_count_test <- 142
print(paste("Number of cases to run:", case_count_test))
print("Look at MaxScore and descide how many candidates to report.")

source("1_import_paths.R")
# source("2_import_load_clean_priors.R")

results <- vector("list", length(file_list))

for (i in seq_along(file_list)) {
  print(paste("importing file", i))
  dt_selected_sub <- readRDS(file_list[i])
  
  print(paste(i, "cleaning"))
  source("2_import_clean.R")

  print(paste(i, "priors"))
  source("3_import_priors.R")
  
  print(paste(i, "formatting"))
  source("4_import_format.R")

  print(paste(i, "storing"))
  results[[i]] <- maindataset
  rm(maindataset)
  gc()
}

maindataset <- bind_rows(results)

saveRDS(maindataset, file = "../data/preprocessed/ACMGuru_singlecase_df_raw_v1_momic_samplecount_180_all.Rds")

maindataset <- readRDS(file = "../data/preprocessed/ACMGuru_singlecase_df_raw_v1_momic_samplecount_180_all.Rds")

maindataset <- maindataset |>
select(all_of(col_names))

# data.table::fwrite(maindataset, file=paste0(path_preprocessed, "ACMGuru_singlecase_df_raw_v1_momic_samplecount_180_all.tsv"), sep = "\t")

maindataset_top <- maindataset |>
  group_by(sample.id) |>
  group_modify(~ {
    if (any(.x$SYMBOL %in% c("MMACHC", "MTHFR"))) {
      slice(.x, 1)
    } else {
      slice_max(.x, ACMG_total_score, n = 1, with_ties = FALSE)
    }
  }) |>
  ungroup() |>
  select(sample.id, ACMG_total_score, SYMBOL, HGVSp, HGVSc, everything())

data.table::fwrite(maindataset_top, file=paste0(path_preprocessed, "ACMGuru_singlecase_df_raw_v1_momic_samplecount_180_top_1.tsv"), sep = "\t")

maindataset_top_10 <- maindataset |>
  group_by(sample.id) |>
  group_modify(~ {
    priority_rows <- dplyr::filter(.x, SYMBOL %in% c("MMACHC", "MTHFR"))
    remaining_rows <- dplyr::filter(.x, !SYMBOL %in% c("MMACHC", "MTHFR"))
    n_remaining <- max(10 - nrow(priority_rows), 0)
    top_remaining <- dplyr::slice_max(remaining_rows, ACMG_total_score, n = n_remaining, with_ties = FALSE)
    dplyr::bind_rows(priority_rows, top_remaining)
  }) |>
  ungroup() |>
  select(sample.id, ACMG_total_score, SYMBOL, HGVSp, HGVSc, everything())

data.table::fwrite(maindataset_top_10, file=paste0(path_preprocessed, "ACMGuru_singlecase_df_raw_v1_momic_samplecount_180_top_10.tsv"), sep = "\t")

maindataset_known_gene <- maindataset |>
  filter(SYMBOL == "MTHFR" | SYMBOL == "MMACHC") |> 
  arrange(sample.id)

data.table::fwrite(maindataset_known_gene, file=paste0(path_preprocessed, "ACMGuru_singlecase_df_raw_v1_momic_samplecount_180_known_gene.tsv"), sep = "\t")


# hold <- dt_selected_sub
# dt_selected_sub <- hold

print("Note that a new 'df' would be sourced in report.Rmd so do not use this variable without checking")

source("5_launch_report.R")



