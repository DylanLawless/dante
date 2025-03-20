library(dplyr)

output_dir <- "../data/"
custom_variable <- ""
df_report_discussion <- readr::read_tsv(file=(paste0(output_dir, "get_discussion/df_report_discussion_", custom_variable, ".tsv")))

df_mim <- data.table::fread( file="../../scicore_mirror/data/acmguru_data_20250107/ACMGuru_singlecase_df_report_v1_momic_samplecount_180_all_chr_rank1_mim.csv", sep = ",")

# run server -----
# install.packages("ollamar")
library(ollamar)
test_connection()
list_models()

# Real world application ----
data <- df_mim |> dplyr::slice(1)

table_1 <- data |> dplyr::select(
"sample.id",
"ACMG_total_score",
"ACMG_count",
"SYMBOL",
"HGVSp",
"HGVSc",
"Consequence",
"IMPACT",
# "genotype",
"gnomAD_AF",
"mim_morbid_description") 

table_1 <- paste(names(table_1), table_1[1, ], sep=": ", collapse=". ")
class(table_1)
intro_text <- "Here is the variant interpretation output from GuRu:"
table_1 <- paste(intro_text, table_1)

# Print the combined text
print(table_1)

gene_discussion <- df_report_discussion %>% dplyr::filter(SYMBOL %in% data$SYMBOL) 
gene_discussion <- paste(names(gene_discussion), gene_discussion[1, ], sep=": ", collapse=". ")

intro_text2 <- "Here is the get_discussion output from GuRu:"
table_2 <- paste(intro_text2, gene_discussion)

# Concatenate 'table_1' and 'gene_discussion' into one string
text_discussion <- paste(table_1, table_2)
text_discussion <- paste(table_2)

# Predefined components -----
text_description <- "
Instructions for Creating the 'Evidence Interpretation Summary' Section in a Clinical Genetics Report:
Title: Begin with the section titled 'Evidence Interpretation Summary'.
Objective: Summarize the biological role of the gene/protein related to the pathogenic variant found in the patient.
Explanation of Clinical Relevance:
Discuss the protein function details and their implications in the disease context.
Highlight the clinical significance of the findings.
Pathogenicity Assessment:
Note that GuRu identified a pathogenic variant for one patient listed in Table 1.
Mention the high total score from GuRu and the complex scoring functions which support the classification of this variant as pathogenic.
Conclude that there is substantial evidence to consider this variant a likely causal genetic determinant of the disease.
Language and Formatting:
Use concise, clinically appropriate language.
Ensure all information provided is directly relevant to the genetic interpretation.
Format the summary to comply with the LaTeX document structure, using appropriate sectioning for clarity and readability.
"

# Complete prompt construction for each loop
prompt <- paste(text_description, text_discussion, sep=" ")

data_source <- "x"
model <- "llama3.2"
model <- "deepseek-r1"

responses <- lapply(1:length(data_source), function(i) {
	start_time <- Sys.time()

	# Generate response
	resp <- generate(model, prompt, stream = TRUE)
	
	# Process and save response
	processed_text <- resp_process(resp, "text")
	writeLines(processed_text, con=paste0("Output_", i, ".txt"))
	
	elapsed_time <- difftime(Sys.time(), start_time, units="secs")
	cat("Model used:", model, "\nTime taken:", elapsed_time, "seconds\n")

	return(processed_text)
})

# Review responses or use them as needed
responses

library(stringr)

processed_text <- responses[1]
parts <- str_split(processed_text, "</think>", 2)
extracted_text <- if(length(parts[[1]]) > 1) trimws(parts[[1]][2]) else ""


# Save the extracted text to file, naming with data source details and model
file_name <- paste0("Output_", data_source, "_", model, ".txt")

writeLines(extracted_text, con=file_name)


# Save the extracted text to file, naming with data source details and model
file_name <- paste0("Output_", data_source[i]$sample_id, "_", model, ".txt")
writeLines(extracted_text, con=file_name)





# [[1]] llama3.2
# [1] "Summary_start\n\n\\section{Biological function summary}\n\nThe patient presents with a pathogenic variant (ACMG.score = 0) in the LILRA1 gene, which encodes for Leukocyte immunoglobulin-like receptor subfamily A member 1 (CD85 antigen-like family member I). The GuRu output suggests that this protein has roles as an inhibitory MHC class I receptor activity [GO:0032396], transmembrane signaling receptor activity [GO:0004888] and acts as receptor for class I MHC antigens.\n\nThis variant, chr19:54596470_T/C, leads to a substitution of Serine (Ser) at position 414 with Proline (Pro), denoted as p.Ser414Pro. The impact of this substitution on the function of the LILRA1 protein is considered moderate.\n\nThe clinical significance of this pathogenic variant in the context of this patient's condition cannot be determined without further clinical evaluation and correlation with phenotype.\n\nSummary_end"

# Model used: deepseek-r1 
# Time taken: 221.3659 seconds
# [1] "<think>\nOkay, so the user is asking me to write a section of a clinical genetics report focusing on variant interpretation in LaTeX format. They provided some specific instructions, so I need to make sure I follow them closely.\n\nFirst, they want it clear and succinct with no unrelated content. That means only relevant information about the variant's interpretation should be included. They emphasized correct reporting of factual evidence, so accuracy is crucial here.\n\nLooking at the data given: a patient detected by GuRu where ACMG.score indicates the variant is pathogenic. The get_discussion output from GuRu provides additional functional annotations for LILRA1. So, I need to include both the genetic assessment and the biological implications.\n\nI should structure this section with bullet points for readability. Each point should address a different aspect—like the variant's effect, its role in the gene, pathway involvement, and implications for diagnosis or management.\n\nWait, the user mentioned using LaTeX, so I'll make sure to format it correctly. The section heading is already provided as \\section{Evidence interpretation summary}, which I'll include at the top. Each bullet point should be concise, starting with the variant's impact, then moving on to its biological role and potential clinical significance.\n\nAlso, I need to ensure that each bullet starts with an asterisk for emphasis but keeps it brief without unnecessary details. The conclusion should tie everything together by reaffirming the variant as pathogenic based on the data.\n\nI think I've covered all the points: genetic assessment, functional annotations, gene pathway involvement, and clinical implications. Now, putting it all together in a structured way that's easy to read and meets the user's requirements.\n</think>\n\n\\section{Evidence interpretation summary}\n* The variant chr19:54596470_T/C (p.Ser414Pro) in LILRA1 is pathogenic based on GuRu's scoring system, with an ACMG.score of MODERATE.  \n* SYMBOL: LILRA1 encodes a protein involved in immunoglobulin-like receptor function, including antigen binding to class I MHC molecules and transmembrane signaling activity.  \n* Pathway involvement includes functions such as inhibitory MHC class I receptor activity and CD85 antigen-like family member I activity.  \n* The variant's missense nature (c.1240T>C) suggests a disruption in protein structure or function critical for maintaining proper gene expression or intercellular signaling.  \n\nSummary_start Summary_end"

# generate("llama3.1", "Tomorrow is a...", output = "df", stream = TRUE)


# About ----
  # Note: You should have at least 8 GB of RAM available to run the 7B models, 16 GB to run the 13B models, and 32 GB to run the 33B models.
  
  #' @article{Lin2025JOSS,
  #' 	author = {Lin, Hause and Safi, Tawab}, 
  #' 	title = {ollamar: An R package for running large language models}, 
  #' 	journal = {Journal of Open Source Software}, 
  #' 	volume = {10}, 
  #' 	number = {105}, 
  #' 	pages = {7211}, 
  #' 	year = {2025},
  #' 	month = jan,
  #' 	volume = {10}, 
  #' 	doi = {10.21105/joss.07211},
  #' 	url = {https://joss.theoj.org/papers/10.21105/joss.07211}
  #' }
  