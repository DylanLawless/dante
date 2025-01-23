library(dplyr)
library(ggplot2)
library(ggrepel)

# library(patchwork)
# library(RColorBrewer)

# # Read subject ID and plot
# singlecase_ID <- "HCY001"

# source: 21b_pca_biplot_1kg_part3_ggplot_vcurrent.R
out_path <- "../../variant_level/data/scicore_PCA/PCA" 
case_out <- "../output/" 

# Example individual case plot for genetic reports ----
df <- readRDS(paste0(out_path, "data_PCA_biplot_swisspedhealth_1kg_singlecase.Rds"))

 
# custom - rename IDs to drop symbol
df <- df %>% mutate(Individual.ID = str_replace_all(Individual.ID, "_", ""))

population_colors <- readRDS(paste0(out_path, "population_colors.Rds"))
legend_order <- readRDS(paste0(out_path, "legend_order.Rds"))
legend_labels <- readRDS(paste0(out_path, "legend_labels.Rds"))
label_count_case <- readRDS(paste0(out_path, "label_count_case.Rds"))
label_count_control <- readRDS(paste0(out_path, "label_count_control.Rds"))

# general plot
p_general <- 
  df |>
  ggplot(aes(x = PC1, y = PC2, color = Population)) +
  geom_point(data = filter(df, Population != "SwissPedHealth"), alpha = 0.5) +
  geom_point(data = filter(df, Population == "SwissPedHealth"), alpha = 0.5) +
  scale_color_manual(values = population_colors, breaks = legend_order, labels = legend_labels) +
  labs(x = "First component", y = "Second component") +
  theme_minimal() +
  guides(alpha = "none", size = "none") +
  labs(title = stringr::str_wrap("Principal component analysis (PCA) of 1000 genomes project, reference genome GRCh38, showing population structure for SwissPedHealth.\n ", 60),
       # caption = paste0("Cases n=", label_count_case, "\nControls n=", label_count_control))
     caption = stringr::str_wrap(
       paste0("Examples:
       CEU Utah residents of European ancestry,
       IBS Iberian populations in Spain; 
MXL Mexican Ancestry in Los Angeles;
PUR Puerto Rican;
PJL Punjabi in Lahore Pakistan;
CLM Colombian in Medellin;
STU Sri Lankan Tamil in the UK;
CHS Han Chinese South.",
              "\nCases n=", label_count_case, "\nControls n=", label_count_control), 60),
     color = "\n\nReference populations") +
  theme(legend.text=element_text(size=6))

# custom patient ID data
p_1kg_swisspedhealth_singlecase <- 
  p_general +
  geom_point(data = filter(df, Individual.ID == singlecase_ID), color = "black", alpha = 0.8) +
  geom_label_repel(
    data = filter(df, Individual.ID == singlecase_ID),
    aes(label = Individual.ID), 
    box.padding = 0.35, point.padding = 0.5, segment.color = 'black',
    size = 3, color = "black", nudge_x = 0.02, nudge_y = 0.02
  ) +
  labs(subtitle = paste0("Subject ID: ", singlecase_ID))

# p_1kg_swisspedhealth_singlecase

ggsave(paste0(case_out, "plot_PCA_biplot_swisspedhealth_1kg_singlecase_", singlecase_ID, ".pdf"), p_1kg_swisspedhealth_singlecase, width = 10, height = 5) 
# ggsave(paste0(case_out, "plot_PCA_biplot_swisspedhealth_1kg_singlecase_", singlecase_ID, ".png"), p_1kg_swisspedhealth_singlecase, width = 10, height = 5) 

