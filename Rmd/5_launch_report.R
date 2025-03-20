
print("Limiting reported hits to case_count_test")
unique_ids <- head(unique_ids, case_count_test)  # limited test set

for(sample_id in unique_ids) {
    
  print(paste("Now running report: ", sample_id))
  print(paste("..............................."))
  current_data <- maindataset[maindataset$sample.id == sample_id, ]
  
  serialized_data <- serialize(current_data, NULL)
  data_checksum <- digest(serialized_data, algo = "md5")
  
  print(paste("Data checksum for sample ID", sample_id, ":", data_checksum))
  
  current_data <- current_data |> unique()
  current_data <- current_data |>
    filter(sample.id == sample.id) |>
    slice(1:10) |>
    as.data.frame()
  
  # current_rank <- current_data$rank |> unique()
  current_rank <- min(current_data$rank)
  print(paste("rank:", current_rank))
  
  report_path <- paste0("../reports/sample_", sample_id , "_report_priority_", current_rank, ".pdf")
  
  render("report.Rmd", output_file = report_path,
         params = list(
           singlecase_ID = sample_id,
           data_checksum = data_checksum,
           pdf_list = pdf_list
         ),
         quiet = TRUE)
  
  log_file <-  paste0("./sample_", sample_id , "_report_priority_", current_rank,  ".log")
  file.remove(log_file)
  print(paste("Done."))
}
