# Make sure that the data match what was used in the analyses.
# If you run this script and it breaks, it indicates the local data
# Do not match the analytic data used to produce the paper

library(tidyverse)

# Load data ----

# NSDUH National 2024 microdata
nsduh <- read_csv("data/nsduh_2024_adults.csv.gz", show_col_types = FALSE)

# State by age from 2023-24 pooled data
state_by_age <- read_csv(
  "data/nsduh_state_by_age_2023_24.csv",
  show_col_types = FALSE
)

# All adults, by state
state_all_adults <- read_csv(
  "data/state_context_all_adults_2023_24.csv",
  show_col_types = FALSE
)

# Checks that underlying data haven't changed ----
stopifnot(
  nrow(nsduh) == 47299,
  n_distinct(nsduh$stratum) == 50,
  n_distinct(paste(nsduh$stratum, nsduh$psu)) == 100,
  setequal(unique(nsduh$age_group), c("18-25", "26-34", "35+")),
  all(na.omit(nsduh$cannabis_dui) %in% c(0, 1)),
  all(na.omit(nsduh$alcohol_dui) %in% c(0, 1)),
  all(nsduh$analysis_weight > 0),
  nrow(state_by_age) == 153,
  all(state_by_age |> count(age_group) |> pull(n) == 51),
  n_distinct(state_by_age$jurisdiction) == 51,
  all(complete.cases(select(
    state_by_age,
    cannabis_dui_estimate, cannabis_dui_ci_lower, cannabis_dui_ci_upper,
    alcohol_dui_estimate, alcohol_dui_ci_lower, alcohol_dui_ci_upper
  ))),
  nrow(state_all_adults) == 51,
  n_distinct(state_all_adults$jurisdiction) == 51,
  all(complete.cases(state_all_adults)),
  sum(state_all_adults$recreational_legal_2024 == 1) == 25,
  sum(state_all_adults$drug_testing_coverage_pct >= 40) == 26
)

message("All input-data checks passed.")
