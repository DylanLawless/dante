output_dir <- "../data/"

df_mim <- data.table::fread( file="../../scicore_mirror/data/acmguru_data_20250107/ACMGuru_singlecase_df_report_v1_momic_samplecount_180_all_chr_rank1_mim.csv", sep = ",")

# get discussion ----

# import ----
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("rtracklayer")
library(rtracklayer)

# uniprot to PDB ----
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")

library(biomaRt)

# listMarts()
# ensembl = useMart("ensembl")
# listDatasets(ensembl)
# ensembl = useDataset('hsapiens_gene_ensembl',mart=ensembl)

# 
# uniprotR ----
# install.packages('UniprotR')
# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("UniprotR", force = TRUE)
library(UniprotR)

# Homo sapiens (Human) [9606] Consists of 204,185 entries

# Data import ----
df_uniprot <- rtracklayer::readGFF("../ref/uniprot/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22.gff")
# df_uniprot <- readGFF("../data/uniprot_HUMAN_Q2TBE0_CWF19L2.gff")
df_uniprot_meta <- read.csv("../ref/uniprot/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22.tab", sep="\t")

# colnames(df)[colnames(df) == 'oldName'] <- 'newName'
colnames(df_uniprot_meta)[colnames(df_uniprot_meta) == 'Entry'] <- 'seqid'

df_uniprot$seqid
df_uniprot_meta$seqid
df_uniprot_meta$Gene.names

# Separate the Gene.names column where names are separated by space
df_uniprot_meta_tidy <- df_uniprot_meta %>% 
  tidyr::separate_rows(Gene.names, sep = " ")

# Rename Gene.names to SYMBOL
colnames(df_uniprot_meta_tidy)[colnames(df_uniprot_meta_tidy) == 'Gene.names'] <- 'SYMBOL'

# get gene list ----
gene_list <- df_mim$SYMBOL |> unique()

df_uniprot_meta_tidy_genes <- 
  df_uniprot_meta_tidy |> 
  filter(SYMBOL %in% gene_list) |>
  filter(Status == "reviewed")

df_uniprot_meta_tidy <- df_uniprot_meta_tidy_genes

rm(df_uniprot, df_uniprot_meta_tidy_genes, df_uniprot_meta)
gc()

Accessions <- df_uniprot_meta_tidy$seqid 

#Get Taxonomy Information 
# TaxaObj <- GetNamesTaxa(Accessions)
# saveRDS(TaxaObj, file=paste(output_dir, "ontology_taxa/TaxaObj.Rds", sep = ""))

#Get Gene ontolgy Information 
# GeneOntologyObj <- GetProteinGOInfo(Accessions)
# saveRDS(GeneOntologyObj, file=paste(output_dir, "ontology_taxa/GeneOntologyObj.Rds", sep = ""))

# GetProteinFunction
# ProteinFunction <- GetProteinFunction(Accessions)
# saveRDS(ProteinFunction, file=paste(output_dir, "ontology_taxa/ProteinFunction.Rds", sep = ""))

# load local copy
TaxaObj <- readRDS(file=paste(output_dir, "ontology_taxa/TaxaObj.Rds", sep = ""))
GeneOntologyObj <- readRDS(file=paste(output_dir, "ontology_taxa/GeneOntologyObj.Rds", sep = ""))
ProteinFunction <- readRDS(file=paste(output_dir, "ontology_taxa/ProteinFunction.Rds", sep = ""))

# get discussion ----
names(GeneOntologyObj)
names(TaxaObj)
names(ProteinFunction)
# GeneOntologyObj[2:10,4]

GeneOntologyObj$seqid  <- rownames(GeneOntologyObj)
TaxaObj$seqid  <- rownames(TaxaObj)
ProteinFunction$seqid  <- rownames(ProteinFunction)

data_discussion <- merge(GeneOntologyObj, TaxaObj)
data_discussion <- merge(data_discussion, ProteinFunction)

# represents top hits 
df_mim_seqid <- df_uniprot_meta_tidy |> dplyr::select(SYMBOL, seqid)
df_report <- merge(df_mim,df_mim_seqid )
data_discussion <- merge(df_report, data_discussion, by = "seqid")

df_report_discussion <-
  # df_report_grouped_grouped_df_max_GO_taxa_funct |>
  data_discussion |>
  dplyr::select(SYMBOL,
                Protein.names,
                "Gene.Ontology..molecular.function.",
                "Function..CC."
  ) |>
  unique()

custom_variable <- ""

readr::write_tsv(df_report_discussion, file=(paste0(output_dir, "get_discussion/df_report_discussion_", custom_variable, ".tsv")))

