library(dplyr)
library(ggplot2)

source("./report_qv3_table.R")

print("Merge with the cohort list of known phenotype")
out_path <- "../../scicore_mirror/data/acmguru_out_20250107/"
data_path <- "../../scicore_mirror/data/acmguru_data_20250107"
getwd()

# data clinical ----
df_clin <- read.csv2(file = paste0(data_path, "/hcy_cohort_patient_data_20240627.csv"))
df_clin$SAMPLE_NAME <- df_clin$sample_id
df_clin <- df_clin %>% dplyr::select(-sample_id)

# convert the dataset -----
df <- dt_selected_sub_rank1
 

df$SAMPLE_NAME <- df$sample.id
df$MAX_AF <- as.numeric(df$MAX_AF)
df <- df %>% mutate(gnomAD_AF = if_else(is.na(gnomAD_AF), 0, gnomAD_AF))
df <- df %>% mutate(MAX_AF = if_else(is.na(MAX_AF), 0, MAX_AF))
df$AF.x <- as.numeric(df$AF.x)

# Merge clinical data ----
df <- merge(df_clin, df, all.y = TRUE)
rm(df_clin)

# drop columns
df <- dplyr::select(df, -c( ACMG_score:ACMG_PP4, Strong_pathogenic_GE:Supporting_pathogenic_GE))

df$ensembl_gene_id <- df$Gene

# Get the best variant ID
df <- df %>% mutate(HGVSpc = coalesce(HGVSp, HGVSc))

df <- df %>% mutate(Symbol_Variant = paste(SYMBOL, HGVSpc, sep = " - ")) 

df <- df %>% dplyr::select("sample_group", "sample_subgroup",
                           "sample.id", 
                           "source", "HGVSpc", "SYMBOL",  
                           "max_score_in_case",
                           # "case_rank_in_cohort",
                           "variant_rank_in_case",
                           "ACMG_total_score",
                           "ACMG_count",
  "genotype", "comp_het_flag", "ensembl_gene_id", "gnomAD_AF", "MAX_AF", "MAX_AF_POPS", everything())

df <- df %>%
  arrange(max_score_in_case ) #%>%
  # arrange(case_rank_in_cohort)





df_mthfr <-  df|> filter(SYMBOL == "MTHFR")
df_mmachc <-  df|> filter(SYMBOL == "MMACHC")
df_tcn2 <-  df|> filter(SYMBOL == "TCN2")
df_mmachc$HGVSpc |> unique()
df_mthfr$HGVSpc |> unique()
df_tcn2$HGVSpc |> unique()

# # expected versus autodetected
# MMACHC	2 variants (55 cases)	 - 2 found variants
# MTHFR (39 cases)	38 variants	- 18 found variants
# TCN2 (1 case) - 1 found variants


## deduplicate ----

# Optimized approach using data.table for better performance
library(data.table)

# Convert df to a data.table
df0 <- df
setDT(df0)

# Adjust the priority settings to handle NAs explicitly
df0[, priority := fifelse(is.na(source), 2, fifelse(source == "metabol", 1, 2))]

# Correctly filtering duplicates with prioritisation
df0 <- df0[order(priority)]  # Order by priority

# Remove duplicates based on all columns except 'source' and 'priority'
df0 <- df0[!duplicated(df0, by = setdiff(names(df0), c("source", "priority")))]

# Optionally, you can remove the 'priority' column if it's no longer needed
df0[, priority := NULL]

df <- df0 %>% as.data.frame()
rm(df0)
names(df)

gc()

# Focus on top scores for now ----
# df_hold <- df
# df <- df %>%
  # filter(ACMG_total_score > 8)

# Biomart ----
# Assuming df is your dataframe
df_unique_genes <- unique(df$ensembl_gene_id)

# Load the biomaRt library
library(biomaRt)

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve MIM Morbid Accession numbers for the unique genes
mim_data <- getBM(attributes = c('ensembl_gene_id', 'mim_morbid_accession', 'mim_gene_accession', 'mim_morbid_description'),
                  filters = 'ensembl_gene_id',
                  values = df_unique_genes,
                  mart = ensembl)

# View the first few entries of the retrieved data
head(mim_data)



# Collapse `mim_data` unique genes
collapsed_mim_data <- mim_data %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    mim_morbid_accession = paste(unique(mim_morbid_accession), collapse = ";"),
    mim_gene_accession = paste(unique(mim_gene_accession), collapse = ";"),
    mim_morbid_description = paste(unique(mim_morbid_description), collapse = ";"),
    .groups = 'drop'
  )


# Merge the MIM data back with the original dataframe
df_mim <- merge(df, collapsed_mim_data, all = TRUE, allow.cartesian = TRUE)

# View the merged dataframe
head(df_mim)
names(df_mim)

df_mim <- df_mim %>% dplyr::select("source", "HGVSpc", "mim_morbid_accession", "mim_morbid_description", "gnomAD_AF", MAX_AF, MAX_AF_POPS, everything())

df_mim <- df_mim %>%
  arrange(max_score_in_case ) #%>%
  # arrange(case_rank_in_cohort)

df_mim$sample_subgroup %>% unique()


# get the diagnosis column ----

library(dplyr)
library(stringr)

df_mim <- df_mim %>%
  mutate(known_diagnosis = if_else(str_detect(sample_subgroup, ":"),
                                   str_extract(sample_subgroup, "^[^:]+"),
                                   sample_subgroup))



df_mim <- df_mim %>%
  mutate(confirmed_diagnosis = if_else(known_diagnosis == SYMBOL, "YES", NA_character_))


df_mim  <- df_mim |> dplyr::select("sample_group", "sample_subgroup",
              "sample.id", 
              "source", "HGVSpc", "SYMBOL",  
              "max_score_in_case",
              # "case_rank_in_cohort",
              "variant_rank_in_case",
              "ACMG_total_score",
              "ACMG_count",
              "genotype", "comp_het_flag", 
              'mim_morbid_accession', 'mim_gene_accession', 'mim_morbid_description',
              "ensembl_gene_id", "gnomAD_AF", "MAX_AF", "MAX_AF_POPS", everything())

df_mim  <- df_mim |> dplyr::select(known_diagnosis, confirmed_diagnosis,  everything())

# if diagnosed drop other variants ----

df_mim <- df_mim %>%
  group_by(sample.id) %>%
  # Filter rows while handling NAs correctly
  filter(if (all(is.na(confirmed_diagnosis))) TRUE else confirmed_diagnosis == "YES") %>%
  ungroup()


# Counting occurrences in known_diagnosis
known_tally <- df_mim %>%
  count(known_diagnosis, sort = TRUE)  # sort = TRUE to order by frequency

# Counting occurrences in confirmed_diagnosis
confirmed_tally <- df_mim %>%
  count(confirmed_diagnosis, sort = TRUE)

# Print the tallies to view them
print(known_tally)
print(confirmed_tally)


library(dplyr)
 
library(dplyr)

# Assuming df_mim already includes the known_diagnosis and confirmed_diagnosis columns
combined_tally <- df_mim %>%
  group_by(known_diagnosis, confirmed_diagnosis) %>%
  summarise(count = n_distinct(sample.id), .groups = 'drop') %>%
  mutate(total = sum(count))  # Calculate total cases across all groups

# Calculate the percentage of total diagnosed cases for each group
combined_tally <- combined_tally %>%
  mutate("cohort%" = (count / total) * 100) %>%
  dplyr::select(-total)  # Remove the total column if it's no longer needed

# Print the combined tally with percentages
print(combined_tally)


25+23+56+19+19+38



df_mim <- df_mim %>%
  arrange(max_score_in_case ) %>%
  arrange(confirmed_diagnosis, known_diagnosis , ACMG_total_score)


df_mim <- df_mim %>% 
  dplyr::select("sample.id", 
                "known_diagnosis", "confirmed_diagnosis",
              # "case_rank_in_cohort", 
              "variant_rank_in_case",
              "ACMG_total_score", "ACMG_count", "source", "SYMBOL", "HGVSp", "HGVSc", 
              "Consequence", "IMPACT", "genotype", "gnomAD_AF", "comp_het_flag", everything())



df_mim <- df_mim |> dplyr::select(sample.id:AF.x)

data.table::fwrite(df_mim, file="../../scicore_mirror/data/acmguru_data_20250107/ACMGuru_singlecase_df_report_v1_momic_samplecount_180_all_chr_rank1_mim.csv", sep = ",")

filter source metabol
undiagnosed

# 
# library(stringr)
# 
# df_mim_cand <- df_mim %>%
#   filter(source == "metabol") %>%
#   filter(MAX_AF < 0.1) %>%
#   filter(str_detect(mim_morbid_description, "METAB"))
# 
# df_mim_cand %>%
#   dplyr::select(SYMBOL, gnomAD_AF, MAX_AF, mim_morbid_description) %>% 
#   unique()
# 
# df_mim_only <- df_mim %>%
#   dplyr::select(SYMBOL, mim_morbid_description) %>% 
#   unique()
# 
# df_mim_cand_sum <-  df_mim_cand %>%
#   group_by(sample_subgroup, SYMBOL, genotype, comp_het_flag, Symbol_Variant, ACMG_total_score, HGVSpc, IMPACT, Consequence, source, mim_morbid_description) %>%
#   summarise(Count = n_distinct(SAMPLE_NAME), .groups = "drop") %>%
#   arrange(SYMBOL, HGVSpc)  # Adjusted to put SYMBOL first for clarity in grouping
# 
# df_mim_cand_sample <-  df_mim_cand %>%
#   group_by(SAMPLE_NAME, genotype, comp_het_flag, SYMBOL, Symbol_Variant, HGVSpc, IMPACT, Consequence, source, mim_morbid_description) %>%
#   summarise(Count = n_distinct(SAMPLE_NAME), .groups = "drop") %>%
#   arrange(SYMBOL, HGVSpc)  # Adjusted to put SYMBOL first for clarity in grouping
# 
# write.csv(df_mim_cand_sum, paste0(out_path,"df_mim_cand_sum", ".csv"))
# write.csv(df_mim_cand_sample, paste0(out_path,"df_mim_cand_sample", ".csv"))
# 
# df_mim_HIBCH <- df_mim %>%
#   filter(SYMBOL == "HIBCH") %>%
#   filter(MAX_AF < 0.1)# %>%
#   # filter(str_detect(mim_morbid_description, "METAB"))
# 
# df_mim$MAX_AF
# 
# # df_mim_cand <- df_mim %>%
# #   # filter(source == "metabol") %>%
# #   filter(gnomAD_AF < 0.001)
# 
# # Undiagnosed ----
#  
# # "MTHFR", "MMACHC:c.394C>T", "MMACHC:c.271dupA", "HCY+_undiagnosed", "HCY+MMA+_undiagnosed"  
# 
# df_mim_undiag <- df_mim %>%
#   filter(sample_subgroup == "HCY+_undiagnosed" | sample_subgroup == "HCY+MMA+_undiagnosed"
# ) 
# 
# # Correctly filtering duplicates with prioritisation
# df_mim_undiag <- df_mim_undiag %>% arrange(desc(ACMG_total_score))
# df_mim_undiag <- df_mim_undiag %>% filter(MAX_AF < 0.1) 
# 
# # remove variants that also occur in diagnosed cases.
# # Subset for the target subgroups
# target_subgroups <- df_mim %>%
#   filter(sample_subgroup %in% c("MTHFR", "MMACHC:c.394C>T", "MMACHC:c.271dupA"))
# 
# # Find HGVSpc values in the target subgroup
# target_hgvs_pc <- target_subgroups$HGVSpc %>% unique()
# 
# # Identify the rows to be filtered out in the undiagnosed subgroups
# df_mim_undiag <- df_mim_undiag %>%
#   filter(!(sample_subgroup %in% c("HCY+_undiagnosed", "HCY+MMA+_undiagnosed") & 
#              HGVSpc %in% target_hgvs_pc))
# 
# df_mim_undiag <- df_mim_undiag %>% filter(source == "metabol")
# 
# write.csv(df_mim_undiag, paste0(out_path,"df_mim_undiag", ".csv"))
# 
# 
# df_mim_undiag_censor <- df_mim_undiag %>%
#   dplyr::select(source:comp_het_flag, ACMG_highest:PUBMED)
# 
# write.csv(df_mim_undiag_censor, paste0(out_path,"df_mim_undiag_censored", ".csv"))
# 
# 
# # %>%
#   # filter(str_detect(mim_morbid_description, "METAB"))
# 
# 
#   
# 
# # Master ACMG score filter ----
# df_hold <- df 
# # df <- df_hold
# df <- df %>% filter(ACMG_total_score >= 14)
# 
# # Get study book data ----
# data_source <- "../../study_book/master_table_cohort/out/df_phase2.Rds"
# df_summary <- readRDS(data_source)
# 
# # cap_source <- paste("Data source:",  "\\newline", sprintf("\\texttt{%s}", gsub("_", "\\_", data_source, fixed = TRUE)))
# # cap_text <- paste("Actual count of raw data. Flag shows whether raw WGS has been processed to gVCF stage. Cohort and data-path come from a shared summary of data expected to be present.", cap_source)
# 
# df_summary$SAMPLE_NAME <- df_summary$Clean_id
# df_summary <- df_summary %>% filter(data_type == "WGS")
# df_summary <- df_summary %>% dplyr::select("Cohort", "data_path", "SAMPLE_NAME" )
# 
# # head(df_summary)
# # head(df_main)
# 
# print("FIX THE SAMPLE NAMES HERE TO MERGE")
# print("FIX THE SAMPLE NAMES HERE TO MERGE")
# 
# df <- merge(df_summary, df)
# 
# # Assuming 'df' has columns SYMBOL, ACMG_total_score, and rownames, and that SAMPLE_NAME identifies unique samples
# df_variant_count <- df %>%
#   group_by(SYMBOL, Symbol_Variant, ACMG_total_score, HGVSpc, genotype, IMPACT, Consequence, source) %>%
#   summarise(Count = n_distinct(SAMPLE_NAME), .groups = "drop") %>%
#   arrange(SYMBOL, HGVSpc)  # Adjusted to put SYMBOL first for clarity in grouping
# 
# # df_variant_count$genotype <- as.character(df_variant_count$genotype)
# # df_variant_count$genotype <- factor(df_variant_count$genotype)
# 
# # df_variant_count$genotype <- factor(df_variant_count$genotype)
# # print(levels(df_variant_count$genotype))  # Optional: check the levels to ensure correctness
# 
# 
# lab_cap <- stringr::str_wrap("[1] metabolic panel: https://panelapp.genomicsengland.co.uk/panels/467/. Likely inborn error of metabolism - targeted testing not possible (Version 5.24), Relevant disorders: Likely inborn error of metabolism - targeted testing not possible, Likely inborn error of metabolism, Inborn errors of metabolism, R98. Panel types: GMS Rare Disease Virtual, Component Of Super Panel, GMS signed-off. Latest signed off version: v5.0 (1 May 2024). [2] IUIS immune disorder panel:  International Union of Immunological Societies (IUIS) Inborn Errors of Immunity Committee (IEI) (https://iuis.org), (Tangye et al., 2022). Cleaned as per https://lawlessgenomics.com/topic/iuis-iei-table-page.", width = 120)
# 
# p1 <- ggplot(df_variant_count, aes(x = Symbol_Variant, y = Count)) +
#   geom_text(aes(label = genotype, color = genotype), vjust = -0.25) +
#   ggrepel::geom_text_repel(aes( y = Count + 5 , label = source, alpha = .4), direction = "y",  nudge_y = 10, angle = 90) +
#   geom_col(aes(fill = ACMG_total_score)) +
#   # geom_text(aes(label = Count), vjust = -0.25) +
#   scale_y_continuous(limits = c(0, max(df_variant_count$Count)+20)) +  # Set y-axis limits
#   scale_colour_gradient(low = "red",high = "purple") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(x = "Gene symbol - Variant", y = "Count of Samples") +
#   guides(alpha = "none") 
# 
# p1_solo <- 
#   p1 +  labs(
#     title = "Count of samples sharing variants with highest pathogenicity score",
#   subtitle = "Genotype (0,1,2)\nCount of samples\nIncl. known gene panels",
#   caption = lab_cap,
#   fill = "ACMG Total Score")
# 
# p1
# p1_solo
# ggsave(filename = paste0(out_path, "/df_variant_count.pdf"), p1_solo, width = 9, height = 7)
# 
# ## gnomad AF ----
# # df$gnomAD_AF, df$MAX_AF, df$MAX_AF_POPS
# df$MAX_AF <- as.numeric(df$MAX_AF)
# df <- df %>% mutate(gnomAD_AF = if_else(is.na(gnomAD_AF), 0, gnomAD_AF))
# df <- df %>% mutate(MAX_AF = if_else(is.na(MAX_AF), 0, MAX_AF))
# 
# # MAX_AF, MAX_AF_PO
# 
# df$AF.x <- as.numeric(df$AF.x)
# 
# p2 <- df %>%
#   ggplot( aes(x = Symbol_Variant, y = -MAX_AF)) +
#   geom_point(aes(size = AF.x, alpha = (1-AF.x))) +
#   # scale_colour_gradient(low = "red",high = "purple") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(x = "Gene symbol - Variant") 
# 
# p2_solo <- p2 + labs(title = "Gnomad frequency, pop max")
# p2_solo
# 
# 
# ggsave(filename = paste0(out_path, "/df_gnomad.pdf"), p2_solo, width = 9, height = 7)
# 
# p2_blank <- 
#   p2 +
#   theme(#axis.line=element_blank(),
#         axis.text.x=element_blank(),
#         # axis.title.x=element_blank()
#         )
# 
# p3 <- p2_blank  + labs(
#   title = "Count of samples sharing variants with highest pathogenicity score",
#   subtitle = "Gnomad and cohort frequency\nGenotype (0,1,2)\nCount of samples\nIncl. known gene panels")
# 
# library(patchwork)
# patch1 <- (p3 /  p1)
# 
# patch1
# 
# ggsave(filename = paste0(out_path, "/df_variant_count_gnomad.pdf"), patch1, width = 9, height = 9)
# 
# 
# ## undiagnosed ----
# 
# df_mim_undiag$mim_morbid_description
# df_mim_undiag$MAX_AF <- as.numeric(df_mim_undiag$MAX_AF)
# df_mim_undiag <- df_mim_undiag %>% mutate(gnomAD_AF = if_else(is.na(gnomAD_AF), 0, gnomAD_AF))
# df_mim_undiag <- df_mim_undiag %>% mutate(MAX_AF = if_else(is.na(MAX_AF), 0, MAX_AF))
# df_mim_undiag$AF.x <- as.numeric(df_mim_undiag$AF.x)
# df_mim_undiag <- df_mim_undiag %>% dplyr::select(AF.y, AF.x, everything())
# 
# # MAX_AF, MAX_AF_PO
# # df$AF.x <- as.numeric(df$AF.x)
# library(dplyr)
# library(stringr)
# 
# df_mim_undiag <- df_mim_undiag %>%
#   mutate(
#     # Convert to sentence case
#     mim_morbid_description = str_to_sentence(mim_morbid_description),
#     # Truncate to 60 characters and create a new column
#     mim_morbid_description_short = str_sub(mim_morbid_description, 1, 20)
#   )
# 
# 
# # Plot
# p4 <- df_mim_undiag %>%
#   ggplot( aes(x = Symbol_Variant, y = -MAX_AF)) +
#   geom_point(aes(size = (AF.x), color = ACMG_total_score)) +
#   scale_colour_gradient(low = "orange",high = "purple") +
#   ggrepel::geom_text_repel(aes( y = -MAX_AF + 0 , label = mim_morbid_description_short, alpha = .4), direction = "y",  nudge_y = 0.05, angle = 90) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(x = "Gene symbol - Variant") 
#   # scale_y_continuous(limits = c(0, max(df_variant_count$Count)+20))
#   
# p4
# 
# # query_terms <- c("cystin", "cobal", "metabolism", "Thiamine", "methylcrotonyl", "dehydrogenase")
# # 
# # # Create a regex pattern that matches any of the terms
# # regex_pattern <- paste(query_terms, collapse = "|")
# # 
# # # Filter rows in df_mim_undiag where mim_morbid_description contains any of the query terms
# # df_mim_undiag_cystin <- df_mim_undiag %>%
# #   filter(str_detect(mim_morbid_description, regex_pattern))
# # 
# # 
# # 
# # df_mim_undiag_cystin %>% dplyr::select(Symbol_Variant, SAMPLE_NAME, mim_morbid_description_short)
# 
# 
# # Define your query terms and create a regex pattern
# query_terms <- c("cystin", "cobal", "metabolism", "Thiamine", "methylcrotonyl", "dehydrogenase")
# regex_pattern <- paste(query_terms, collapse = "|")
# 
# # Add a new column to the data frame to indicate if the description matches any query terms
# df_mim_undiag <- df_mim_undiag %>%
#   mutate(match_highlight = ifelse(str_detect(mim_morbid_description, regex_pattern), 1, 0))
# 
# 
# # Plotting
# p4 <- ggplot(df_mim_undiag, aes(x = Symbol_Variant, y = -MAX_AF)) +
#   geom_point(aes(size = AF.x, color = ACMG_total_score)) +
#   scale_colour_gradient(low = "orange", high = "purple") +
#   ggrepel::geom_text_repel(aes(y = -MAX_AF + 0, label = mim_morbid_description_short, alpha = match_highlight * 0.7 + 0.3),
#                            direction = "y", nudge_y = 0.05, angle = 90,
#                            size = 3,  # Adjust text size
#                            max.overlaps = Inf) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(x = "Gene symbol - Variant")
# 
# # Print the plot
# p4
# 
# 
# lab_cap <- stringr::str_wrap("[1] metabolic panel: https://panelapp.genomicsengland.co.uk/panels/467/. Likely inborn error of metabolism - targeted testing not possible (Version 5.24), Relevant disorders: Likely inborn error of metabolism - targeted testing not possible, Likely inborn error of metabolism, Inborn errors of metabolism, R98. Panel types: GMS Rare Disease Virtual, Component Of Super Panel, GMS signed-off. Latest signed off version: v5.0 (1 May 2024). All phenotypes: MTHFR, MMACHC:c.394C>T, MMACHC:c.271dupA, HCY+_undiagnosed, HCY+MMA+_undiagnosed", width = 120)
# 
# lab_title <- "Undiagnosed samples with highest pathogenicity variants."
# lab_subtitle <- "All known metabolic gene panel\nGenotype (0,1,2)\nCount of samples\nIncl. known gene panels\nNote fixed PP3 scores!"
# 
# p4 <- p4 + labs(title = lab_title, subtitle = lab_subtitle, caption = lab_cap)
#   
# ggsave(filename = paste0(out_path, "/df_mim_undiag.pdf"), p4, width = 9, height = 9)
# 
# # TCN2 ----
# # In the remaining 60 samples, we identified bi-allelic disease- causing variants for 22 individuals in genes other than MMUT: ACSF3 (17 individuals) (Fig. 1h), TCN2 (3 individuals), SUCLA2 (1 individual) and MMAB (1 individual) according to the ACMG classification (Supple- mentary Table 1). 
# 
# df_hold %>% filter(SYMBOL == "ACSF3")
# x <- df_hold %>% na.omit(source)
# x <- df_hold %>% filter(source == "iuis" | source == "metabol")
# 
# x <- df_hold %>% filter(ACMG_total_score > 6)
# 
# 
# gene_list <- c("ACSF3", "TCN2", "SUCLA2", "MMAB", "MMACHC")
# df_gene_list <- df_hold %>%
#   filter(SYMBOL %in% gene_list ) %>% 
#   group_by(SYMBOL, Symbol_Variant, 
#            # ACMG_total_score, 
#            HGVSpc, genotype, IMPACT, Consequence) %>%
#   summarise(Count = n_distinct(SAMPLE_NAME), .groups = "drop") 
# 
# df_gene_list_long <- df_hold %>%
#   filter(SYMBOL %in% gene_list ) %>% 
#   group_by(SAMPLE_NAME, SYMBOL, Symbol_Variant, ACMG_total_score, HGVSpc, genotype, IMPACT, Consequence) %>%
#   summarise(Count = n_distinct(SAMPLE_NAME), .groups = "drop") 
# 
# 
# 
# # 
# # # Biomart test
# # # if (!requireNamespace("biomaRt", quietly = TRUE)) {
# # # install.packages("BiocManager")
# # # BiocManager::install("biomaRt")
# # # }
# # 
# # # vargen 
# # # https://github.com/MCorentin/vargen
# # # 
# # # install.packages("devtools")
# # library(devtools)
# # # install_github(repo = "MCorentin/vargen", dependencies = TRUE)
# # library(vargen)
# # 
# # # vargen_install(install_dir = "./vargen_data", gtex_version = "v8", verbose = T)
# # 
# # gene_mart <- connect_to_gene_ensembl()
# # View(list_omim_accessions(gene_mart))
# # 
# # # You can search using list of keywords as well:
# # list_omim_accessions(gene_mart, c("alzheimer", "neurodegeneration"))
# # 
# # 
# # gene_mart <- connect_to_gene_ensembl()
# # View(list_omim_accessions(gene_mart))
# # 
# # # You can search using list of keywords as well:
# # list_omim_accessions(gene_mart, c("alzheimer", "neurodegeneration"))
# # 
# # #>      mim_morbid_accession    mim_morbid_description
# # #> [1]  606889                  ALZHEIMER DISEASE 4;;AD4;;ALZHEIMER DISEASE, FAMILIAL, 4
# # #> [2]  618476                  BRAIN ABNORMALITIES, NEURODEGENERATION, AND DYSOSTEOSCLEROSIS;BANDDOS
# # #> [3]  615643                  NEURODEGENERATION WITH BRAIN IRON ACCUMULATION 6; NBIA6
# # #> [4]  618276                  NEURODEGENERATION, CHILDHOOD-ONSET, WITH CEREBELLAR ATROPHY; CONDCA
# # #> [5]  618278                  FIBROSIS, NEURODEGENERATION, AND CEREBRAL ANGIOMATOSIS; FINCA
# # #> ...
# # 
# # 
# # 
# # 
# # # Biomart test
# # # if (!requireNamespace("biomaRt", quietly = TRUE)) {
# # # install.packages("BiocManager")
# # # BiocManager::install("biomaRt")
# # # }
# # library(biomaRt)
# # # 
# # # # Connect to Ensembl
# # # ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# # # 
# # # # Retrieve gene IDs (and additional attributes as needed)
# # # genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position', 'strand'), mart = ensembl)
# # # 
# # # # View the first few entries of the retrieved data
# # # head(genes)
# # # 
# # # # Retrieve gene IDs and MIM Morbid Accession
# # # genes_with_mim <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'mim_morbid_accession', 'mim_morbid_description'),
# # #                         mart = ensembl)
# # # 
# # # # View the first few entries of the retrieved data
# # # head(genes_with_mim)
# # # 
# # # 
# # # 
# # # df$Ensembl
# # # df$Gene
# # # 
# # # 
# # # searchDatasets(mart = ensembl, pattern = "hsapiens")
# # # 
# # # 
# # # 
# # # # Load the biomaRt package
# # # library(biomaRt)
# # # 
# # # # Connect to the Ensembl BioMart database for human genes
# # # ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# # # 
# # # # Specify the gene name for which you want to retrieve data
# # # gene_name <- "TP53"
# # # 
# # # # Define the attributes you want to retrieve
# # # # attributes <- c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "mim_morbid_description")
# # # 
# # # attributes <- c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position")
# # # 
# # # # Retrieve data for the specified gene
# # # gene_data <- getBM(attributes = attributes,
# # #                    filters = "hgnc_symbol",
# # #                    values = gene_name,
# # #                    mart = ensembl)
# # # 
# # # # Display the retrieved data
# # # print(gene_data)
# # # 
# # # chromosome_name <- "21"
# # # 
# # # # Retrieve data for all genes on chromosome 1
# # # genes_chromosome_21 <- getBM(attributes = attributes,
# # #                              filters = "chromosome_name",
# # #                              values = chromosome_name,
# # #                              mart = ensembl)
# # # 
# # # # Display the retrieved data
# # # print(genes_chromosome_1)
# # # 
# # 
# # 
