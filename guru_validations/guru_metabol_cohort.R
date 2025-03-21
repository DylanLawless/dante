library(readxl)
library(dplyr)
library(ggplot2);theme_set(theme_bw())
library(tidyr)
library(stringr)
library(ggrepel)

df <- readxl::read_excel("../data/preprocessed/ACMGuru_singlecase_df_raw_v1_momic_samplecount_180_known_gene.xlsx")
df$gnomAD_AF <-  as.numeric(df$gnomAD_AF)

df$sample.number <- as.numeric(sub("HCY", "", df$sample.id))

max(df$gnomAD_AF)

missing_numbers <- setdiff(1:150, df$sample.number)
df <- bind_rows(df, data.frame(sample.number = missing_numbers, opinion = 0, genotype_name = "NULL"))
df$gnomAD_AF[is.na(df$gnomAD_AF)] <- 0

df |> 
  filter(gnomAD_AF < 0.05) |>
  ggplot(aes(x = sample.number, y = opinion, fill = genotype_name)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~ SYMBOL, ncol = 1, scales = "free_y") +
  labs(title  = "QV with gnomad_AF < 0.05 (max 0.01)")

max(df$gnomAD_AF)

p1 <- df |> 
  ggplot(aes(x = sample.number, y = opinion, fill = genotype_name)) +
  geom_vline(xintercept = 40.5, linetype = "dotted", color = "red") +
  geom_vline(xintercept = 100.5, linetype = "dotted", color = "red") +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~ SYMBOL, ncol = 1, scales = "free_y") +
  labs(title  = "QV with with gnomad_AF < 0.05 (max gnomAD 0.3)",
       subtitle = "Shows the 7 excluded samples in 1:100",
       caption = str_wrap(paste("Stacks indicate multiple variants in a subject, such as compound heterozygous. Opinion indicates ease of interpretation, e.g. 0 = benign, 1-3 = heterozygous missense to homozygous LoF.")), 40)

p1 
ggsave(file = "../data/guru_validations/p_guru_metabol_cohort.pdf", p1, height = 6, width = 7)


# validation ----
df2 <- read.table("../data/validation_phase2/phase2_deepvariant_filtered_commented_20240607_known.csv", sep = ",", header = TRUE)


names(df)
names(df2)

# colnames(df2)[colnames(df2) == 'oldName'] <- 'newName'
colnames(df2)[colnames(df2) == 'VCFRowID'] <- 'rownames'
colnames(df2)[colnames(df2) == 'sample_id'] <- 'sample.id'
df2$sample.number <- as.numeric(sub("HCY", "", df2$sample.id))
missing_numbers <- setdiff(1:150, df2$sample.number)
df2 <- bind_rows(df2, data.frame(sample.number = missing_numbers, opinion = 1, genotype_name = "NULL"))
df$gnomAD_AF[is.na(df$gnomAD_AF)] <- 0

df2$source <- "deepvar"
df$source <- "gatk"

df_select <- df |> filter(gnomAD_AF < 0.05) |> select(source, sample.number, rownames) # HGVSp, HGVSc, 
df2_select <- df2 |> select(source, sample.number, rownames) # HGVSp, HGVSc, 

head(df_select$rownames)
head(df2_select$rownames)
df_select$rownames <- gsub("[: /]", "_", df_select$rownames)

df2_select$opinion <- 1

head(df_select) 
head(df2_select)

df_gatk <- df_select %>% 
  filter(!is.na(rownames)) %>% 
  distinct(sample.number, rownames) %>% 
  mutate(gatk = TRUE)

df_deepvar <- df2_select %>% 
  filter(!is.na(rownames)) %>% 
  distinct(sample.number, rownames) %>% 
  mutate(deepvar = TRUE)

df_all <- full_join(df_gatk, df_deepvar, by = c("sample.number", "rownames")) %>%
  mutate(category = case_when(
    !is.na(gatk) & !is.na(deepvar) ~ "validated",
    is.na(gatk) & !is.na(deepvar) ~ "missing",
    !is.na(gatk) & is.na(deepvar) ~ "novel"
  ))

missing_labels <- df_all %>%
  filter(sample.number < 101) %>%
  filter(category == "missing") %>%
  group_by(sample.number) %>%
  summarise(missing_count = n(), .groups = "drop") %>%
  mutate(ypos = missing_count / 2)

p2 <- ggplot(df_all, aes(x = sample.number)) +
  geom_vline(xintercept = 40.5, linetype = "dotted", color = "red") +
  geom_vline(xintercept = 100.5, linetype = "dotted", color = "red") +
  geom_bar(position = "stack", aes( fill = category), color = "black") +
  geom_text_repel(
    data = missing_labels,
    aes(x = (sample.number), y = ypos+.5, label = sample.number),
    nudge_y = 2,
    segment.color = "grey10"
  ) +
  labs(x = "sample.number", y = "count of variants")

p2
ggsave(file = "../data/guru_validations/p_guru_metabol_cohort_validation.pdf", p2, height = 6, width = 8)


print("Is the first missing variant present in the gVCF?")

# system('gzcat ../data/validation_phase2/chr1_bcftools_gatk_norm_decom_maf.recode_vep_chr1_11010153_12010152.vcf.gz | grep "11801196" | grep "1/1" | grep "ENSP00000365669.3:p.Gln147Pro"')


# system('gzcat ../data/validation_phase2/chr1_bcftools_gatk_norm_decom_maf.recode_vep_dbnsfp_pick.vcf.gz | grep "45507544"')







print("Systematic check for missing variant present in the gVCF?")

results <- character(0)
for(rn in unique(df_missing$rownames)) {
  num <- strsplit(rn, "_")[[1]][2]
  # cmd <- paste0('gzcat ../data/validation_phase2/chr1_bcftools_gatk_norm_decom_maf.recode_vep_chr1_11010153_12010152.vcf.gz | grep "', num, '"')
  
  cmd <- paste0('grep ../data/validation_phase2/chr1_bcftools_gatk_norm_decom_maf.recode_vep_dbnsfp_pick.vcf "', num, '" | head -n1')
  
  output <- system(cmd, intern = TRUE)
  results <- c(results, paste0(rn, ": ", length(output), " lines"))
}
caption_text <- paste(results, collapse = "\n")
print(caption_text)

cap <- paste("A known variant is missing for 16 samples\n",
             "This consistes of", length(unique(df_missing$rownames)),
             "unique variants.\n", 
             "Variants present in gVCF:\n",
             caption_text)


cap
 
for term in 11801196 11793907 11790916 11794385 11803155 11794029 45507544 45508329; do count=$(zcat /project/data/shared_all/release_dna_snv_indel_v2/annotation/chr1_bcftools_gatk_norm_decom_maf.recode_vep_dbnsfp_pick.vcf.gz | grep -m 1 "$term" | wc -l); echo "$term: $count"; done

capvar <- paste("11801196: 1
11793907: 1
11790916: 1
11794385: 1
11803155: 1
11794029: 1
45507544: 1
45508329: 1")


cap <- paste("A known variant is missing for 16 samples\n",
             "This consistes of", length(unique(df_missing$rownames)),
             "unique variants.\n", 
             "Variants present in gVCF:\n",
             capvar)

cap

p3 <- p2 + labs(caption = cap)
p3

ggsave(file = "../data/guru_validations/p_guru_metabol_cohort_validation.pdf", p3, height = 6, width = 8)



