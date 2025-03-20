# dt_selected_sub <- file_list |>
#   lapply(readRDS) |>
#   bind_rows()

results <- vector("list", length(file_list))

for (i in seq_along(file_list)) {
  print(paste("importing file", i))
  dt_selected_sub <- readRDS(file_list[i])
  source("import_clean.R")
  source("import_format.R")
  source("import_priors.R")
  results[[i]] <- maindataset
  rm(maindataset)
  gc()
}

maindataset <- bind_rows(results)

