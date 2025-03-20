
colnames(dt_selected_sub)[colnames(dt_selected_sub) == 'sample'] <- 'sample.id'

# only ~5% of positions have carriers
dt_selected_sub <- dt_selected_sub |> filter(genotype > 0) 
# dt_selected_sub <- dt_selected_sub |> filter(gnomAD_AF < 0.5)

dt_selected_sub <- dt_selected_sub |>
  filter(is.na(gnomAD_AF) | gnomAD_AF < 0.5)  # painful lesson learned here chasing these NA 

# dt_selected_sub <- dt_selected_sub |> filter(ACMG_total_score > 9)

dt_selected_sub <- dt_selected_sub %>%
  filter(str_starts(sample.id, "HCY"))

dt_selected_sub <- dt_selected_sub %>%
  mutate(sample.id_meta = str_extract(sample.id, "_NGS.*"),
         sample.id = str_remove(sample.id, "_NGS.*"))

names(dt_selected_sub)

dt_selected_sub <- dt_selected_sub |> 
  select("sample.id", 
         "ACMG_PVS1":"ACMG_PP4",
         "rownames",
         "ACMG_total_score",
         "ACMG_count",
         "SYMBOL",
         "HGVSp", 
         "HGVSc",
         "Consequence",
         "IMPACT",
         "genotype", 
         "gnomAD_AF",
         "comp_het_flag"  
         # "ACMG_score", "ACMG_highest", "chr", "cohort_pheno","frequency_in_1", Inheritance", "CLIN_SIG","SymbolCount"
  )

data.table::fwrite(dt_selected_sub, file=paste0(path_clean_unprocessed, "ACMGuru_singlecase_df_raw_v1_momic_samplecount_180_chunk_", i, ".tsv"), sep = "\t")
# 
dt_selected_sub_impact_high <- dt_selected_sub |> filter(IMPACT == "HIGH") |> select("sample.id", 
       "rownames",
       "ACMG_total_score",
       "SYMBOL",
       "HGVSp", 
       "HGVSc",
       "Consequence",
       "IMPACT",
       "genotype", 
       "gnomAD_AF",
       "comp_het_flag") 

data.table::fwrite(dt_selected_sub_impact_high, file=paste0(path_clean_unprocessed_impact, "ACMGuru_singlecase_df_raw_v1_momic_samplecount_180_high_impact_chunk_", i, ".tsv"), sep = "\t")

rm(dt_selected_sub_impact_high)
