print("Adding priors to dt_selected_sub")
dt_selected_sub <- dt_selected_sub |>
  mutate(ACMG_total_score = if_else(SYMBOL %in% c("TCN2", "MMACHC", "MTHFR"), ACMG_total_score + prior_score_weight, ACMG_total_score))

