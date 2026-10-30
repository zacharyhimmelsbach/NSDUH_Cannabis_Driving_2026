# Reproduce the results, tables, and figures in manuscript.

# Load libraries
library(tidyverse)
library(survey)
library(sandwich)
library(stargazer)

#########################################
# Setup ----
#########################################

## Tell survey package to error if any cells only have 1 unit
options(survey.lonely.psu = "fail")

## Parameters for age ranges and substance use frequency categories
age_levels <- c("18-25", "26-34", "35+")
frequency_levels <- c(
  "<1/week", "1 to <4/week", "4 to <6/week", "Approximately 6-7/week"
)

#####################################
## Helper Functions ----
#####################################

### Format p-values
format_p <- function(x) {
  if_else(x < 0.001, "less than .001", str_remove(sprintf("%.3f", x), "^0"))
}

### Stop stargazer from printing
quiet_stargazer <- function(...) {
  invisible(capture.output(stargazer(...)))
}

### Calculate confidence intervals for proportion of reporters (table 1)
proportion_ci <- function(estimate, standard_error, degrees_freedom) {
  critical <- qt(0.975, degrees_freedom)
  logit_se <- standard_error / (estimate * (1 - estimate))
  c(
    lower = plogis(qlogis(estimate) - critical * logit_se),
    upper = plogis(qlogis(estimate) + critical * logit_se)
  )
}

### Wrapper to get HC3 standard errors
hc3_results <- function(model, population) {
  beta <- coef(model)
  covariance <- vcovHC(model, type = "HC3")
  standard_error <- sqrt(diag(covariance))
  degrees_freedom <- df.residual(model)
  critical <- qt(0.975, degrees_freedom)

  tibble(
    population = population,
    term = names(beta),
    estimate_percentage_points = 100 * unname(beta),
    standard_error_percentage_points = 100 * standard_error,
    ci_lower_percentage_points = 100 * (beta - critical * standard_error),
    ci_upper_percentage_points = 100 * (beta + critical * standard_error),
    p_value = 2 * pt(-abs(beta / standard_error), degrees_freedom),
    jurisdictions = nobs(model)
  )
}

#########################################
# Load and clean data ----
#########################################

# National microdata
nsduh <- read_csv("data/nsduh_2024_adults.csv.gz", show_col_types = FALSE) |>
  mutate(
    age_group = factor(age_group, levels = age_levels),
    sex = factor(sex, levels = c("Female", "Male")),
    race_ethnicity = factor(
      race_ethnicity,
      levels = c(
        "Non-Hispanic White", "Non-Hispanic Black", "Hispanic",
        "Non-Hispanic Asian", "Non-Hispanic other/multiple races"
      )
    ),
    poverty_group = factor(
      poverty_group,
      levels = c(">=200% FPL", "100-199% FPL", "<100% FPL")
    ),
    metro_status = factor(
      metro_status,
      levels = c("Large metro", "Small metro", "Nonmetro")
    ),
    frequency_group = factor(
      case_when(
        past_year_cannabis == 1 & cannabis_days_past_year < 52 ~ "<1/week",
        past_year_cannabis == 1 & cannabis_days_past_year < 208 ~ "1 to <4/week",
        past_year_cannabis == 1 & cannabis_days_past_year < 300 ~ "4 to <6/week",
        past_year_cannabis == 1 & cannabis_days_past_year <= 365 ~
          "Approximately 6-7/week",
        TRUE ~ NA_character_
      ),
      levels = frequency_levels
    ),
    complete_pair = !is.na(cannabis_dui) & !is.na(alcohol_dui),
    dui_difference = if_else(complete_pair, cannabis_dui - alcohol_dui, NA_real_),
    model_complete_case = past_year_cannabis == 1 &
      !is.na(cannabis_dui) & !is.na(frequency_group) & !is.na(sex) &
      !is.na(race_ethnicity) & !is.na(poverty_group) & !is.na(metro_status)
  )

# State-by-age 2023-24 pooled data
state_by_age <- read_csv(
  "data/nsduh_state_by_age_2023_24.csv",
  show_col_types = FALSE
) |>
  mutate(age_group = factor(age_group, levels = age_levels))

# State-level aggregates (all adults), 2023-24 pooled data
state_all_adults <- read_csv(
  "data/state_context_all_adults_2023_24.csv",
  show_col_types = FALSE
)

#########################################
# Table 1: national cannabis and alcohol DUI prevalence by age ----
#########################################

national_design <- svydesign(
  ids = ~psu,
  strata = ~stratum,
  weights = ~analysis_weight,
  data = nsduh,
  nest = TRUE
)
design_df <- degf(national_design)

table_1 <- map_dfr(age_levels, function(age) {
  paired_design <- subset(
    national_design,
    age_group == age & complete_pair
  )
  estimates <- svymean(
    ~cannabis_dui + alcohol_dui + dui_difference,
    paired_design,
    na.rm = TRUE
  )
  estimate <- coef(estimates)
  standard_error <- SE(estimates)
  cannabis_ci <- proportion_ci(
    estimate[["cannabis_dui"]],
    standard_error[["cannabis_dui"]],
    design_df
  )
  alcohol_ci <- proportion_ci(
    estimate[["alcohol_dui"]],
    standard_error[["alcohol_dui"]],
    design_df
  )
  difference_critical <- qt(0.975, design_df)

  outcome_design <- subset(national_design, age_group == age)
  cannabis_outcome_mean <- coef(
    svymean(~cannabis_dui, outcome_design, na.rm = TRUE)
  )[[1]]
  alcohol_outcome_mean <- coef(
    svymean(~alcohol_dui, outcome_design, na.rm = TRUE)
  )[[1]]
  age_population <- nsduh |>
    filter(age_group == age) |>
    summarise(total = sum(analysis_weight)) |>
    pull(total)

  tibble(
    age_group = age,
    cannabis_dui_pct = 100 * estimate[["cannabis_dui"]],
    cannabis_ci_lower = 100 * cannabis_ci[["lower"]],
    cannabis_ci_upper = 100 * cannabis_ci[["upper"]],
    alcohol_dui_pct = 100 * estimate[["alcohol_dui"]],
    alcohol_ci_lower = 100 * alcohol_ci[["lower"]],
    alcohol_ci_upper = 100 * alcohol_ci[["upper"]],
    difference_pct_points = 100 * estimate[["dui_difference"]],
    difference_ci_lower = 100 * (
      estimate[["dui_difference"]] -
        difference_critical * standard_error[["dui_difference"]]
    ),
    difference_ci_upper = 100 * (
      estimate[["dui_difference"]] +
        difference_critical * standard_error[["dui_difference"]]
    ),
    relative_cannabis_vs_alcohol_pct = 100 * (
      estimate[["cannabis_dui"]] / estimate[["alcohol_dui"]] - 1
    ),
    p_value = 2 * pt(
      -abs(estimate[["dui_difference"]] / standard_error[["dui_difference"]]),
      design_df
    ),
    estimated_cannabis_dui_persons = age_population * cannabis_outcome_mean,
    estimated_alcohol_dui_persons = age_population * alcohol_outcome_mean,
    unweighted_complete_pairs = nsduh |>
      filter(age_group == age, complete_pair) |>
      nrow()
  )
})

write_csv(table_1, "results/table_1_national_prevalence.csv")

# Make formatted html version (this is what appears in the paper)
table_1_display <- table_1 |>
  transmute(
    `Age group` = age_group,
    `Cannabis, % (95% CI)` = sprintf(
      "%.1f (%.1f-%.1f)", cannabis_dui_pct, cannabis_ci_lower, cannabis_ci_upper
    ),
    `Alcohol, % (95% CI)` = sprintf(
      "%.1f (%.1f-%.1f)", alcohol_dui_pct, alcohol_ci_lower, alcohol_ci_upper
    ),
    `Difference, pp (95% CI)` = sprintf(
      "%.1f (%.1f to %.1f)",
      difference_pct_points, difference_ci_lower, difference_ci_upper
    ),
    `Cannabis vs alcohol, %` = sprintf("%.1f", relative_cannabis_vs_alcohol_pct),
    `P value` = format_p(p_value),
    `Cannabis DUI, persons` = scales::comma(round(estimated_cannabis_dui_persons, -3)),
    `Alcohol DUI, persons` = scales::comma(round(estimated_alcohol_dui_persons, -3))
  )

quiet_stargazer(
  as.data.frame(table_1_display),
  type = "html",
  summary = FALSE,
  rownames = FALSE,
  title = "Table 1. National Prevalence of Self-Reported Driving Under the Influence by Age Group",
  out = "results/table_1_national_prevalence.html"
)

#########################################
# Table 2: survey-weighted cannabis-DUI model (among cannabis users) ----
#########################################

user_design <- subset(national_design, model_complete_case)
user_model <- svyglm(
  cannabis_dui ~ age_group + frequency_group + sex + race_ethnicity +
    poverty_group + metro_status,
  design = user_design,
  family = quasibinomial()
)

user_beta <- coef(user_model)
user_se <- SE(user_model)
user_df <- degf(user_design)
user_critical <- qt(0.975, user_df)
user_p <- 2 * pt(-abs(user_beta / user_se), user_df)
user_ci <- cbind(
  exp(user_beta - user_critical * user_se),
  exp(user_beta + user_critical * user_se)
)

table_2 <- tibble(
  term = names(user_beta),
  adjusted_odds_ratio = exp(user_beta),
  ci_lower = user_ci[, 1],
  ci_upper = user_ci[, 2],
  p_value = user_p,
  design_degrees_freedom = user_df,
  unweighted_n = nrow(user_design$variables)
)

write_csv(table_2, "results/table_2_cannabis_dui_predictors.csv")

# Formatted version of table (this appears in the paper)
quiet_stargazer(
  user_model,
  type = "html",
  title = "Table 2. Predictors of Self-Reported Cannabis DUI Among Past-Year Cannabis Users",
  column.labels = "Adjusted odds ratio",
  covariate.labels = c(
    "Age 26-34 vs 18-25", "Age 35+ vs 18-25",
    "1 to fewer than 4 days/week vs fewer than 1 day/week",
    "4 to fewer than 6 days/week vs fewer than 1 day/week",
    "Approximately 6-7 days/week vs fewer than 1 day/week",
    "Male vs female",
    "Non-Hispanic Black vs non-Hispanic White",
    "Hispanic vs non-Hispanic White",
    "Non-Hispanic Asian vs non-Hispanic White",
    "Other/multiple races vs non-Hispanic White",
    "100%-199% FPL vs at least 200% FPL",
    "Under 100% FPL vs at least 200% FPL",
    "Small metro vs large metro", "Nonmetro vs large metro"
  ),
  coef = list(exp(user_beta)),
  p = list(user_p),
  ci = TRUE,
  ci.custom = list(user_ci),
  report = "vc*s",
  omit.stat = c("ll", "aic"),
  digits = 2,
  notes = "Exponentiated coefficients are adjusted odds ratios; parentheses contain 95% confidence intervals. Survey design df = 50.",
  notes.append = FALSE,
  out = "results/table_2_cannabis_dui_predictors.html"
)

# Marginal probabilities quoted in the Results section.
user_model_data <- user_design$variables
adjusted_probabilities <- map_dfr(frequency_levels, function(frequency) {
  new_data <- user_model_data |>
    mutate(
      age_group = factor("18-25", levels = age_levels),
      frequency_group = factor(frequency, levels = frequency_levels)
    )
  tibble(
    age_group = "18-25",
    frequency_group = frequency,
    adjusted_probability_pct = 100 * weighted.mean(
      predict(user_model, newdata = new_data, type = "response"),
      new_data$analysis_weight
    )
  )
})
write_csv(
  adjusted_probabilities,
  "results/table_2_adjusted_probabilities_18_25.csv"
)

#########################################
# Table 3: variation in state cannabis-DUI rates (by age) ----
#########################################

table_3 <- state_by_age |>
  group_by(age_group) |>
  summarise(
    jurisdictions = n(),
    mean_pct = 100 * mean(cannabis_dui_estimate),
    median_pct = 100 * median(cannabis_dui_estimate),
    iqr_lower_pct = 100 * quantile(cannabis_dui_estimate, 0.25),
    iqr_upper_pct = 100 * quantile(cannabis_dui_estimate, 0.75),
    minimum_pct = 100 * min(cannabis_dui_estimate),
    minimum_jurisdiction = jurisdiction[which.min(cannabis_dui_estimate)],
    maximum_pct = 100 * max(cannabis_dui_estimate),
    maximum_jurisdiction = jurisdiction[which.max(cannabis_dui_estimate)],
    .groups = "drop"
  )

write_csv(table_3, "results/table_3_state_distribution.csv")

# Formatted version (this appears in the manuscript)
table_3_display <- table_3 |>
  transmute(
    `Age group` = as.character(age_group),
    `Mean, %` = sprintf("%.1f", mean_pct),
    `Median, %` = sprintf("%.1f", median_pct),
    `IQR, %` = sprintf("%.1f-%.1f", iqr_lower_pct, iqr_upper_pct),
    `Range, % (jurisdiction)` = sprintf(
      "%.1f (%s) to %.1f (%s)",
      minimum_pct, minimum_jurisdiction, maximum_pct, maximum_jurisdiction
    )
  )

quiet_stargazer(
  as.data.frame(table_3_display),
  type = "html",
  summary = FALSE,
  rownames = FALSE,
  title = "Table 3. State-Level Distribution of Self-Reported Cannabis-Impaired Driving",
  out = "results/table_3_state_distribution.html"
)

#########################################
# Appendix Table 1: state cannabis-DUI prevalence by age ----
#########################################

appendix_table_1 <- state_by_age |>
  transmute(
    jurisdiction,
    age_group,
    estimate = sprintf(
      "%.1f (%.1f-%.1f)",
      100 * cannabis_dui_estimate,
      100 * cannabis_dui_ci_lower,
      100 * cannabis_dui_ci_upper
    )
  ) |>
  pivot_wider(names_from = age_group, values_from = estimate) |>
  arrange(jurisdiction) |>
  rename(
    State = jurisdiction,
    `18-25 years, % (95% CI)` = `18-25`,
    `26-34 years, % (95% CI)` = `26-34`,
    `35+ years, % (95% CI)` = `35+`
  )

write_csv(
  appendix_table_1,
  "results/appendix_table_1_state_cannabis_dui.csv"
)

quiet_stargazer(
  as.data.frame(appendix_table_1),
  type = "html",
  summary = FALSE,
  rownames = FALSE,
  title = "Appendix Table 1. State-Level Estimates of Self-Reported Cannabis-Impaired Driving by Age Group, 2023-24",
  out = "results/appendix_table_1_state_cannabis_dui.html"
)

#########################################
# Table 4: state daily-use and recreational-legalization models
#########################################

state_18_25 <- state_by_age |>
  filter(age_group == "18-25") |>
  mutate(daily_use_per_5pp = daily_cannabis_use_estimate / 0.05)

state_model_18_25 <- lm(
  cannabis_dui_estimate ~ recreational_legal_2024 + daily_use_per_5pp,
  data = state_18_25
)
state_model_all_adults <- lm(
  cannabis_dui_prevalence ~ recreational_legal_2024 + daily_use_per_5pp,
  data = state_all_adults
)

table_4 <- bind_rows(
  hc3_results(state_model_18_25, "Ages 18-25"),
  hc3_results(state_model_all_adults, "Adults aged 18 or older")
)
write_csv(table_4, "results/table_4_state_predictors.csv")

state_models <- list(state_model_18_25, state_model_all_adults)
state_model_tables <- list(
  table_4 |> filter(population == "Ages 18-25"),
  table_4 |> filter(population == "Adults aged 18 or older")
)
state_coefs <- map(state_models, ~ 100 * coef(.x))
state_p <- map(state_model_tables, ~ .x$p_value)
state_ci <- map(
  state_model_tables,
  ~ as.matrix(select(.x, ci_lower_percentage_points, ci_upper_percentage_points))
)

quiet_stargazer(
  state_model_18_25,
  state_model_all_adults,
  type = "html",
  title = "Table 4. State-Level Predictors of Cannabis DUI Prevalence",
  column.labels = c("Ages 18-25", "Adults 18+"),
  dep.var.labels.include = FALSE,
  covariate.labels = c(
    "Recreational cannabis legal by end of 2024",
    "Daily/near-daily cannabis use prevalence (per 5 pp)"
  ),
  coef = state_coefs,
  p = state_p,
  ci = TRUE,
  ci.custom = state_ci,
  report = "vc*s",
  omit.stat = c("ser", "f"),
  digits = 2,
  add.lines = list(c("HC3 robust inference", "Yes", "Yes")),
  notes = "Coefficients and 95% confidence intervals are percentage-point differences across 51 jurisdictions.",
  notes.append = FALSE,
  out = "results/table_4_state_predictors.html"
)

#########################################
# FARS/toxicology Spearman validity check reported in the manuscript ----
#########################################

fars_analysis <- state_all_adults |>
  filter(drug_testing_coverage_pct >= 40)

spearman_test <- cor.test(
  fars_analysis$cannabis_dui_prevalence,
  fars_analysis$cannabinoid_positive_pct_among_tested_drivers,
  method = "spearman",
  exact = FALSE
)

fars_results <- tibble(
  method = "Spearman",
  correlation = unname(spearman_test$estimate),
  ci_lower = NA_real_,
  ci_upper = NA_real_,
  p_value = spearman_test$p.value,
  states = nrow(fars_analysis),
  minimum_testing_coverage_pct = 40
)
write_csv(fars_results, "results/fars_external_validity.csv")

# Recalculate Spearman's rho across state toxicology-reporting cutoffs.
fars_spearman_by_cutoff <- map_dfr(0:60, function(cutoff) {
  cutoff_data <- state_all_adults |>
    filter(drug_testing_coverage_pct >= cutoff)
  cutoff_test <- cor.test(
    cutoff_data$cannabis_dui_prevalence,
    cutoff_data$cannabinoid_positive_pct_among_tested_drivers,
    method = "spearman",
    exact = FALSE
  )

  tibble(
    minimum_testing_coverage_pct = cutoff,
    spearman_rho = unname(cutoff_test$estimate),
    p_value = cutoff_test$p.value,
    states = nrow(cutoff_data)
  )
})
write_csv(
  fars_spearman_by_cutoff,
  "results/fars_spearman_by_reporting_cutoff.csv"
)

#########################################
# Figures: state estimates and FARS sensitivity analysis
#########################################

# Helper function to produce same estimates for both substances
make_state_figure <- function(substance, color, title) {
  estimate_name <- paste0(substance, "_dui_estimate")
  lower_name <- paste0(substance, "_dui_ci_lower")
  upper_name <- paste0(substance, "_dui_ci_upper")

  plot_data <- state_by_age |>
    transmute(
      age_group,
      jurisdiction,
      estimate = .data[[estimate_name]],
      lower = .data[[lower_name]],
      upper = .data[[upper_name]],
      state_key = paste(age_group, jurisdiction, sep = "___")
    ) |>
    arrange(age_group, estimate) |>
    mutate(state_key = factor(state_key, levels = unique(state_key)))

  ggplot(plot_data, aes(x = 100 * estimate, y = state_key)) +
    geom_segment(
      aes(x = 100 * lower, xend = 100 * upper, yend = state_key),
      color = "grey75",
      linewidth = 0.35
    ) +
    geom_point(color = color, size = 1.6) +
    facet_wrap(~age_group, scales = "free_y", nrow = 1) +
    scale_y_discrete(labels = ~str_remove(.x, ".*___")) +
    scale_x_continuous(
      breaks = scales::breaks_width(10),
      labels = scales::label_percent(scale = 1, accuracy = 1),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    labs(
      title = title,
      subtitle = "Points show prevalence; horizontal lines show published 95% confidence intervals",
      x = "Prevalence",
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 5.5),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

figure_1 <- make_state_figure(
  "cannabis",
  "#D95F02",
  "State Estimates of Self-Reported Cannabis-Impaired Driving, 2023-24"
)
appendix_figure_1 <- make_state_figure(
  "alcohol",
  "#1F78B4",
  "State Estimates of Self-Reported Alcohol-Impaired Driving, 2023-24"
)
appendix_figure_2 <- ggplot(
  fars_spearman_by_cutoff,
  aes(x = minimum_testing_coverage_pct, y = spearman_rho)
) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.4) +
  geom_vline(
    xintercept = 40,
    color = "grey45",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_step(color = "#5E3C99", linewidth = 0.8) +
  geom_point(
    data = fars_spearman_by_cutoff |>
      filter(minimum_testing_coverage_pct %% 10 == 0),
    color = "#5E3C99",
    size = 2
  ) +
  scale_x_continuous(breaks = seq(0, 60, 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Spearman Correlation Across State Toxicology Reporting Cutoffs",
    subtitle = str_wrap(
      paste(
        "Correlation between self-reported cannabis DUI prevalence and",
        "cannabinoid positivity among tested fatal-crash drivers"
      ),
      width = 80
    ),
    x = "Minimum state drug-testing coverage required (%)",
    y = expression("Spearman " * rho),
    caption = "Dashed line indicates the 40% cutoff used in the primary analysis."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# Save each as png and pdf
ggsave(
  "results/figure_1_state_cannabis_dui.png",
  figure_1,
  width = 10,
  height = 10,
  dpi = 300
)
ggsave(
  "results/figure_1_state_cannabis_dui.pdf",
  figure_1,
  width = 10,
  height = 10
)
ggsave(
  "results/appendix_figure_1_state_alcohol_dui.png",
  appendix_figure_1,
  width = 10,
  height = 10,
  dpi = 300
)
ggsave(
  "results/appendix_figure_1_state_alcohol_dui.pdf",
  appendix_figure_1,
  width = 10,
  height = 10
)
ggsave(
  "results/appendix_figure_2_spearman_by_reporting_cutoff.png",
  appendix_figure_2,
  width = 7,
  height = 5,
  dpi = 300
)
ggsave(
  "results/appendix_figure_2_spearman_by_reporting_cutoff.pdf",
  appendix_figure_2,
  width = 7,
  height = 5
)

#########################################
# Combined table of main results discussed in paper ----
#########################################

main_results <- bind_rows(
  table_1 |>
    transmute(
      result = paste("Cannabis minus alcohol DUI, age", age_group),
      estimate = difference_pct_points,
      ci_lower = difference_ci_lower,
      ci_upper = difference_ci_upper,
      p_value,
      unit = "percentage points"
    ),
  table_2 |>
    filter(term %in% c(
      "age_group26-34", "age_group35+",
      "frequency_groupApproximately 6-7/week", "sexMale",
      "poverty_group100-199% FPL", "poverty_group<100% FPL"
    )) |>
    transmute(
      result = paste("Cannabis-user model:", term),
      estimate = adjusted_odds_ratio,
      ci_lower,
      ci_upper,
      p_value,
      unit = "adjusted odds ratio"
    ),
  table_4 |>
    filter(term != "(Intercept)") |>
    transmute(
      result = paste(population, term),
      estimate = estimate_percentage_points,
      ci_lower = ci_lower_percentage_points,
      ci_upper = ci_upper_percentage_points,
      p_value,
      unit = "percentage points"
    ),
  fars_results |>
    transmute(
      result = paste("FARS external validity:", method),
      estimate = correlation,
      ci_lower,
      ci_upper,
      p_value,
      unit = "correlation"
    )
)
write_csv(main_results, "results/main_results.csv")

message("Finished. Tables and figures are in results folder.")
