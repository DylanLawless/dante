# install.packages("ollamar")
library(ollamar)

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

# example use ----
test_connection()  # test connection to Ollama server
# if you see "Ollama local server not running or wrong server," Ollama app/server isn't running

# download a model
# pull("llama3.1")  # download a model (equivalent bash code: ollama run llama3.1)

# generate a response/text based on a prompt; returns an httr2 response by default
resp <- generate("llama3.2", "tell me a 5-word story")
resp

#' interpret httr2 response object
#' <httr2_response>
#' POST http://127.0.0.1:11434/api/generate  # endpoint
#' Status: 200 OK  # if successful, status code should be 200 OK
#' Content-Type: application/json
#' Body: In memory (414 bytes)

# get just the text from the response object
resp_process(resp, "text")
# get the text as a tibble dataframe
resp_process(resp, "df")

# alternatively, specify the output type when calling the function initially
txt <- generate("llama3.1", "tell me a 5-word story", output = "text")

# list available models (models you've pulled/downloaded)
list_models()


# test ----
resp <- generate("deepseek-r1", "tell me a 5-word story")
resp
resp_process(resp, "text")


# Real world application ----

# Predefined components
description_text <- "You are writing the variant interpretation step of a clinical genetics report in LaTeX format. It must be clear and succinct and therefore have no unrelated content so that the reader can understand. The is important to report the provided factual evidence correctly. The final output should start and end with the machine-readable string Summary_start Summary_end so that we can format automatically. This final result will be in LaTeX format with the section heading \\section{Evidence interpretation summary}."

# Table data and gene discussion information could be sourced from external data handling
variant_info <- "Here is a genetic variant for 1 patient in table 1 as detected by GuRu. Based on the ACMG.score and the numerous complex scoring functions of GuRu we say that this variant is pathogenic."

table_1 <- "Table 1: sample.id pathway_id SYMBOL Variant.GRCh38 HGVSp HGVSc Consequence IMPACT genotype Inheritance gnomAD_AF ACMG.score 261 22 LILRA1 chr19:54596470_T/C ENSP00000251372.3:p.Ser414Pro ENST00000251372.8:c.1240T>C missense_variant MODERATE 1 NA 0 8"

gene_discussion <- "Here is the get_discussion output from GuRu for this gene: SYMBOL Protein.names Gene.Ontology..molecular.function. Function..CC. LILRA1 Leukocyte immunoglobulin-like receptor subfamily A member 1 (CD85 antigen-like family member I) (Leukocyte immunoglobulin-like receptor 6) (LIR-6) (CD antigen CD85i) antigen binding [GO:0003823]; inhibitory MHC class I receptor activity [GO:0032396]; transmembrane signaling receptor activity [GO:0004888] FUNCTION: May act as receptor for class I MHC antigens."

# Complete prompt construction for each loop
prompt <- paste(description_text, variant_info, table_1, gene_discussion, sep=" ")

data_source <- "x"
model <- "llama3.2"
model <- "deepseek-r1"

# Simulating multiple entries handling
# Assuming 'generate' function and 'resp_process' are defined to handle the generation and processing of responses
responses <- lapply(1:length(data_source), function(i) {
	# Start timing
	start_time <- Sys.time()
	
	# Assuming data_source contains different gene data for multiple patients
	# Modify prompt with actual data from data_source[i]
	modified_prompt <- prompt  # Here you would modify the prompt according to data_source[i]
	
	# Generate response
	resp <- generate(model, modified_prompt, stream = TRUE)
	
	# Process and save response
	processed_text <- resp_process(resp, "text")
	writeLines(processed_text, con=paste0("Output_", i, ".txt"))
	
	# Calculate elapsed time
	elapsed_time <- difftime(Sys.time(), start_time, units="secs")
	
	# Print model and time taken
	cat("Model used:", model, "\nTime taken:", elapsed_time, "seconds\n")
	

	return(processed_text)
})

# Review responses or use them as needed
responses



library(stringr)


processed_text <- responses[1]

# Split the text at "</think>" and take everything after it
parts <- str_split(processed_text, "</think>", 2)
extracted_text <- if(length(parts[[1]]) > 1) trimws(parts[[1]][2]) else ""



# Save the extracted text to file, naming with data source details and model
file_name <- paste0("Output_", data_source[i]$sample_id, "_", model, ".txt")
writeLines(extracted_text, con=file_name)

# Save the extracted text to file, naming with data source details and model
file_name <- paste0("Output_", data_source[i]$sample_id, "_", model, ".txt")
writeLines(extracted_text, con=file_name)


<think>
</think>



# [[1]] llama3.2
# [1] "Summary_start\n\n\\section{Biological function summary}\n\nThe patient presents with a pathogenic variant (ACMG.score = 0) in the LILRA1 gene, which encodes for Leukocyte immunoglobulin-like receptor subfamily A member 1 (CD85 antigen-like family member I). The GuRu output suggests that this protein has roles as an inhibitory MHC class I receptor activity [GO:0032396], transmembrane signaling receptor activity [GO:0004888] and acts as receptor for class I MHC antigens.\n\nThis variant, chr19:54596470_T/C, leads to a substitution of Serine (Ser) at position 414 with Proline (Pro), denoted as p.Ser414Pro. The impact of this substitution on the function of the LILRA1 protein is considered moderate.\n\nThe clinical significance of this pathogenic variant in the context of this patient's condition cannot be determined without further clinical evaluation and correlation with phenotype.\n\nSummary_end"

# Model used: deepseek-r1 
# Time taken: 221.3659 seconds
# [1] "<think>\nOkay, so the user is asking me to write a section of a clinical genetics report focusing on variant interpretation in LaTeX format. They provided some specific instructions, so I need to make sure I follow them closely.\n\nFirst, they want it clear and succinct with no unrelated content. That means only relevant information about the variant's interpretation should be included. They emphasized correct reporting of factual evidence, so accuracy is crucial here.\n\nLooking at the data given: a patient detected by GuRu where ACMG.score indicates the variant is pathogenic. The get_discussion output from GuRu provides additional functional annotations for LILRA1. So, I need to include both the genetic assessment and the biological implications.\n\nI should structure this section with bullet points for readability. Each point should address a different aspect—like the variant's effect, its role in the gene, pathway involvement, and implications for diagnosis or management.\n\nWait, the user mentioned using LaTeX, so I'll make sure to format it correctly. The section heading is already provided as \\section{Evidence interpretation summary}, which I'll include at the top. Each bullet point should be concise, starting with the variant's impact, then moving on to its biological role and potential clinical significance.\n\nAlso, I need to ensure that each bullet starts with an asterisk for emphasis but keeps it brief without unnecessary details. The conclusion should tie everything together by reaffirming the variant as pathogenic based on the data.\n\nI think I've covered all the points: genetic assessment, functional annotations, gene pathway involvement, and clinical implications. Now, putting it all together in a structured way that's easy to read and meets the user's requirements.\n</think>\n\n\\section{Evidence interpretation summary}\n* The variant chr19:54596470_T/C (p.Ser414Pro) in LILRA1 is pathogenic based on GuRu's scoring system, with an ACMG.score of MODERATE.  \n* SYMBOL: LILRA1 encodes a protein involved in immunoglobulin-like receptor function, including antigen binding to class I MHC molecules and transmembrane signaling activity.  \n* Pathway involvement includes functions such as inhibitory MHC class I receptor activity and CD85 antigen-like family member I activity.  \n* The variant's missense nature (c.1240T>C) suggests a disruption in protein structure or function critical for maintaining proper gene expression or intercellular signaling.  \n\nSummary_start Summary_end"

# generate("llama3.1", "Tomorrow is a...", output = "df", stream = TRUE)

