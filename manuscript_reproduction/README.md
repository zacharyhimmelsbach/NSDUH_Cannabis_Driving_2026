############################################################ 

# Replication package for

# "**Cannabis-Impaired Driving in the United States: Age- and State-specific Evidence, 2023-2024**"

############################################################ 

The data and code in this folder reproduce the numerical results, tables, and figures that appear in the manuscript.

## Run

From this folder:

``` sh
check_data.R # Checks that input data have not been modified
run_analysis.R # Reproduces results, tables, and figures
```

The scripts require the following R packages: `tidyverse`, `survey`, `sandwich`, and `stargazer`.

## Contents

-   `run_analysis.R`: the single main script. It creates all results, tables, and figures.

-   `check_data.R`: input shape, completeness, and survey-design checks.

-   `data/`: the three analytic files used by the main script. See the readme inside that folder for details.

-   `results/`: CSV and HTML tables; PNG/PDF figures.
