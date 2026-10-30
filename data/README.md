# Analytic data

-   `nsduh_2024_adults.csv.gz`: adult records from the 2024 NSDUH public-use file, reduced to the survey-design, outcome, cannabis-use, and model-covariate fields needed by `run_analysis.R`. DUI outcomes use the NSDUH-provided `DRVINMARJ2` and `DRVINALCO2` recodes (which just collapse all the different categories of answers indicating they did not self-report DUI into a single "No" code).
-   `nsduh_state_by_age_2023_24.csv`: 153 pooled DAS records (51 jurisdictions by 3 age groups) with cannabis- and alcohol-DUI estimates and confidence intervals.
-   `state_context_all_adults_2023_24.csv`: 51 all-adult jurisdiction records with cannabis/alcohol DUI prevalence, daily cannabis use, 2024 recreational-legalization status, and aggregated FARS toxicology-testing metrics.
