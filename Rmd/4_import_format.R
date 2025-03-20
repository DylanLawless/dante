
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

# dt_selected_sub$sample.id_long <- dt_selected_sub$sample.id
# dt_selected_sub$sample.id_long %>% unique()

# Splitting the string based on "NGS" and creating a new column
# dt_selected_sub <- dt_selected_sub %>%
#   mutate(sample.id_meta = str_extract(sample.id, "_NGS.*"),
#          sample.id = str_remove(sample.id, "_NGS.*"))

# dt_selected_sub$sample.id %>% unique()
# dt_selected_sub$sample.id_meta %>% unique()

maindataset <- dt_selected_sub %>%
  mutate(
    HGVSc = str_replace(HGVSc, ":", ": "),
    HGVSp = str_replace(HGVSp, ":", ": "),
    # Consequence = str_replace_all(Consequence, "_", "- ")
    Consequence = str_replace_all(Consequence, "_", " "),
    Consequence = str_replace_all(Consequence, "&", " and ") # required to escape the & in latex tables
  )

maindataset <- maindataset %>%
  mutate(genotype_name = case_when(
    genotype == 1 ~ "Heterozygous",
    genotype == 2 ~ "Homozygous",
    TRUE ~ as.character(genotype)
  ))

names(maindataset)

maindataset <- maindataset %>%
  rowwise() %>%
  mutate(ACMG_criteria = paste0(
    c(
      if(ACMG_PVS1 > 0) "PVS1" else NULL,
      if(ACMG_PS1 > 0) "PS1" else NULL,
      if(ACMG_PS2 > 0) "PS2" else NULL,
      if(ACMG_PS3 > 0) "PS3" else NULL,
      if(ACMG_PS4 > 0) "PS4" else NULL,
      if(ACMG_PS5 > 0) "PS5" else NULL,
      if(ACMG_PM1 > 0) "PM1" else NULL,
      if(ACMG_PM2 > 0) "PM2" else NULL,
      if(ACMG_PM3 > 0) "PM3" else NULL,
      if(ACMG_PM4 > 0) "PM4" else NULL,
      if(ACMG_PM5 > 0) "PM5" else NULL,
      if(ACMG_PM6 > 0) "PM6" else NULL,
      if(ACMG_PM7 > 0) "PM7" else NULL,
      if(ACMG_PP1 > 0) "PP1" else NULL,
      if(ACMG_PP2 > 0) "PP2" else NULL,
      if(ACMG_PP3 > 0) "PP3" else NULL,
      if(ACMG_PP4 > 0) "PP4" else NULL
    ),
    collapse = ", "
  )) %>%
  ungroup()

# mutate(across(c(HGVSc, HGVSp), ~str_wrap(.x, 10))) %>%
# mutate(across(everything(), ~gsub("\n", "<br> ", .x))) %>% # Replace \n with LaTeX \\ for line breaks
maindataset <- maindataset %>%
  group_by(sample.id) %>%
  mutate(MaxScore = max(ACMG_total_score)) %>%
  
  # Filter to keep top 5
  slice_max(order_by = ACMG_total_score, n = 5, with_ties = TRUE) %>%
  
  # Filter to keep only the rows where the score equals the maximum score
  # filter(ACMG_total_score == MaxScore) %>% 
  
  # Ungroup to perform operations not confined by the group
  ungroup() %>%
  # Count the number of unique SYMBOLs per sample.id
  add_count(sample.id, SYMBOL, name = "SymbolCount") %>%
  # Filter to keep only those groups where there's more than one SYMBOL
  group_by(sample.id) %>%
  # filter(any(SymbolCount > 1)) %>% 
  ungroup()

# replace chars for latex ----
replace_chars_for_latex <- function(x) {
  stringr::str_replace_all(x, c("&" = "and", "_" = " ", "\\+" = "$+$"))
}

maindataset <- maindataset %>%
  mutate_if(is.character, replace_chars_for_latex)

maindataset <- maindataset %>% arrange(rank)

# colnames(df)[colnames(df) == 'oldName'] <- 'newName'
# colnames(maindataset)[colnames(maindataset) == 'ACMG_total_score'] <- 'ACMG total score'
# colnames(maindataset)[colnames(maindataset) == 'ACMG_count'] <- 'ACMG count'
# colnames(maindataset)[colnames(maindataset) == 'comp_het_flag'] <- 'comp het flag'

# Deduplicate to get unique sample.ids for report generation
unique_ids <- unique(maindataset$sample.id) 
